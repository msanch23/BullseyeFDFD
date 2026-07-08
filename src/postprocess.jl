function show_sim(ϵ, grid;
                  N_step = 0,
                  flux_mon_top     = ((grid.Z - grid.PMLz_top) + grid.cbg_top) / 2,
                  flux_mon_right   = (grid.total_cbg_radius + (grid.R - grid.PMLρ)) / 2,
                  flux_mon_bot     = (grid.PMLz_bot + grid.cbg_bot) / 2,
                  farfield_mon_top = flux_mon_top,
                  x_limits=nothing, y_limits=nothing)
    xl = isnothing(x_limits) ? (0, grid.R) : x_limits
    yl = isnothing(y_limits) ? (0, grid.Z) : y_limits

    p = heatmap(grid.ρp, grid.zp, abs.(sqrt.(ϵ))',
                aspect_ratio = :equal,
                xlabel = "ρ (μm)",
                ylabel = "z (μm)",
                title = "Bullseye Cavity",
                legend = :bottomright,
                xlims = xl,
                ylims = yl,
                grid = false,
                colormap = cgrad(:grays, rev=true))

    if N_step > 0
        vline!(p, grid.ρp[1:N_step:end], color=:gray, lw=0.5, alpha=0.5, label=false)
        hline!(p, grid.zp[1:N_step:end], color=:gray, lw=0.5, alpha=0.5, label=false)
    end

    if grid.PMLρ > 0
        vspan!(p, [grid.R - grid.PMLρ, grid.R], color=:red, alpha=0.1, label="PML")
        vline!(p, [grid.R - grid.PMLρ], color=:red, lw=1.0, label=false)
    end
    if grid.PMLz_bot > 0
        hspan!(p, [0, grid.PMLz_bot], color=:red, alpha=0.1, label=false)
        hline!(p, [grid.PMLz_bot], color=:red, lw=1.0, label=false)
    end
    if grid.PMLz_top > 0
        hspan!(p, [grid.Z - grid.PMLz_top, grid.Z], color=:red, alpha=0.1, label=false)
        hline!(p, [grid.Z - grid.PMLz_top], color=:red, lw=1.0, label=false)
    end

    has_flux = !isnothing(flux_mon_top) || !isnothing(flux_mon_right) || !isnothing(flux_mon_bot)
    if has_flux
        idx_top   = isnothing(flux_mon_top)   ? nothing : argmin(abs.(grid.zd .- flux_mon_top))
        idx_right = isnothing(flux_mon_right) ? nothing : argmin(abs.(grid.ρd .- flux_mon_right))
        idx_bot   = isnothing(flux_mon_bot)   ? nothing : argmin(abs.(grid.zd .- flux_mon_bot))

        z_top_val = isnothing(idx_top)   ? nothing : grid.zd[idx_top]
        r_val     = isnothing(idx_right) ? nothing : grid.ρd[idx_right]
        z_bot_val = isnothing(idx_bot)   ? nothing : grid.zd[idx_bot]

        r_end = isnothing(r_val) ? grid.R - grid.PMLρ : r_val

        flbl = "Flux Monitor"
        if !isnothing(z_top_val)
            plot!(p, [0, r_end], [z_top_val, z_top_val],
                  color=:green, lw=2.5, linestyle=:dash, label=flbl); flbl = false
        end
        if !isnothing(z_bot_val)
            plot!(p, [0, r_end], [z_bot_val, z_bot_val],
                  color=:green, lw=2.5, linestyle=:dash, label=flbl); flbl = false
        end
        if !isnothing(r_val)
            z_start = isnothing(z_bot_val) ? 0.0 : z_bot_val
            z_end   = isnothing(z_top_val) ? grid.Z : z_top_val
            plot!(p, [r_val, r_val], [z_start, z_end],
                  color=:green, lw=2.5, linestyle=:dash, label=flbl)
        end
    end

    if !isnothing(farfield_mon_top)
        ff_idx_top = argmin(abs.(grid.zd .- farfield_mon_top))
        ff_z_top = grid.zd[ff_idx_top]
        ff_r_end = grid.R - grid.PMLρ
        plot!(p, [0, ff_r_end], [ff_z_top, ff_z_top],
              color=:blue, lw=2.5, linestyle=:dot, label="Far-field Monitor")
    end

    display(p)
    return p
end

cell_dV(grid, i, j) = π * (grid.ρp[i+1]^2 - grid.ρp[i]^2) * grid.Δzp[j]

function _energy_density(mode, grid, inv_ϵr)
    N = grid.N
    ε_ρρ = reshape(1.0 ./ diag(inv_ϵr)[1:N],     grid.Nρ, grid.Nz)
    ε_ϕϕ = reshape(1.0 ./ diag(inv_ϵr)[N+1:2N],  grid.Nρ, grid.Nz)
    ε_zz = reshape(1.0 ./ diag(inv_ϵr)[2N+1:3N], grid.Nρ, grid.Nz)
    return real.(ε_ρρ) .* abs2.(mode.Eρ) .+
           real.(ε_ϕϕ) .* abs2.(mode.Eϕ) .+
           real.(ε_zz) .* abs2.(mode.Ez)
end

function confinement_factor(W, grid)
    int_eye = 0.0; int_cbg = 0.0
    for j in grid.idx_cbg_bot:grid.idx_cbg_top, i in 1:grid.Nρ
        grid.ρd[i] <= grid.total_cbg_radius || continue
        val = W[i, j] * cell_dV(grid, i, j)
        int_cbg += val
        grid.ρd[i] <= grid.eye_radius && (int_eye += val)
    end
    return int_eye / int_cbg
end

function analyze_mode(mode, grid, ϵ, inv_ϵr)
    W = _energy_density(mode, grid, inv_ϵr)

    idx_cbg_r = findlast(grid.ρd .<= grid.total_cbg_radius)

    total_energy = 0.0
    for j in grid.idx_cbg_bot:grid.idx_cbg_top, i in 1:idx_cbg_r
        total_energy += W[i, j] * cell_dV(grid, i, j)
    end
    max_W = maximum(W[1:idx_cbg_r, grid.idx_cbg_bot:grid.idx_cbg_top])

    CF       = confinement_factor(W, grid)
    n_cavity = sqrt(maximum(real.(ϵ[1, grid.idx_cbg_bot:grid.idx_cbg_top])))
    V_eff    = (total_energy / max_W) / (mode.λ / n_cavity)^3
    Fp       = (3.0 / (4.0*π^2)) * (mode.Q / V_eff)

    E_squared = abs2.(mode.Eρ) .+ abs2.(mode.Eϕ) .+ abs2.(mode.Ez)
    normE     = sqrt.(E_squared)

    return (; λ=mode.λ, Q=mode.Q, V_eff, Fp, CF,
              raw_E=mode.raw_E, mode.Eρ, mode.Eϕ, mode.Ez,
              E_squared, normE, W, n_cavity, k_sq=mode.k_sq)
end

function source_power(grid, driven)
    η0 = 376.73
    ρi = driven.source.ρ_idx
    zi = driven.source.z_idx
    dV   = 2π * grid.ρd[ρi] * grid.Δρd[ρi] * grid.Δzd[zi]
    flat = (zi - 1) * grid.Nρ + ρi
    s = conj(driven.J[flat])            * driven.raw_E[flat]            +
        conj(driven.J[flat + grid.N])   * driven.raw_E[flat + grid.N]   +
        conj(driven.J[flat + 2*grid.N]) * driven.raw_E[flat + 2*grid.N]
    return 0.5 * real(s) * dV / η0
end

function analyze_driven(driven, grid)
    power = source_power(grid, driven)
    return (; driven..., power, Purcell=NaN)
end

function report_modes(modes)
    isempty(modes) && return
    has_coll = hasproperty(first(modes), :η)
    w  = has_coll ? 70 : 49
    println("\nEigenmode Analysis Report"); println(repeat("=", w))
    if has_coll
        @printf("%-4s | %-7s | %-6s | %-5s | %-6s | %-6s | %-7s | %-6s\n",
                "mode","λ (nm)","Q","Fp","η (%)","Gauss.","V(λ/n)³","CF")
    else
        @printf("%-4s | %-7s | %-6s | %-5s | %-7s | %-6s\n",
                "mode","λ (nm)","Q","Fp","V(λ/n)³","CF")
    end
    println(repeat("-", w))
    for (i, m) in enumerate(modes)
        if has_coll
            @printf("%-4d | %-7.2f | %-6.1f | %-5.1f | %-6.2f | %-6.4f | %-7.4f | %-6.4f\n",
                    i, m.λ*1000, m.Q, m.Fp, 100*m.η, m.gaussicity, m.V_eff, m.CF)
        else
            @printf("%-4d | %-7.2f | %-6.1f | %-5.1f | %-7.4f | %-6.4f\n",
                    i, m.λ*1000, m.Q, m.Fp, m.V_eff, m.CF)
        end
    end
    println(repeat("=", w)); println()
end

function report_driven(driven)
    if isnan(driven.Purcell)
        @printf("Purcell factor = N/A (purcell=false) @ %.2f nm\n", driven.λ*1000)
    else
        @printf("Purcell factor = %.3f @ %.2f nm\n", driven.Purcell, driven.λ*1000)
    end
    @printf("%s polarized dipole @ (ρ,z): (%d, %d)\n",
            string(driven.source.pol), driven.source.ρ_idx, driven.source.z_idx)
end

function calc_H_fields(grid, mode, Ce; μ=1)
    k0     = 2π / mode.λ
    raw_H  = (1im / 376.73) .* (Ce * mode.raw_E) ./ (μ * k0)
    Hρ = reshape(raw_H[1:grid.N],            grid.Nρ, grid.Nz)
    Hϕ = reshape(raw_H[grid.N+1:2grid.N],    grid.Nρ, grid.Nz)
    Hz = reshape(raw_H[2grid.N+1:3grid.N],   grid.Nρ, grid.Nz)
    return raw_H, Hρ, Hϕ, Hz
end

function calc_farfield(grid, mode, Ce, top_farfield=nothing; m=1, resolution=1080)
    k0 = 2π / mode.λ
    Z0 = 376.73
    _, Hρ_full, Hϕ_full, _ = calc_H_fields(grid, mode, Ce)

    ff_z  = something(top_farfield, ((grid.Z - grid.PMLz_top) + grid.cbg_top) / 2)
    idx_z = argmin(abs.(grid.zp .- ff_z))

    avg(v)  = 0.5 .* (v[1:end-1] .+ v[2:end])
    Eρ_line = mode.Eρ[1:end-1, idx_z]
    Eϕ_line = avg(mode.Eϕ[:, idx_z])
    Hρ_line = avg(0.5 .* (Hρ_full[:, idx_z-1] .+ Hρ_full[:, idx_z]))
    Hϕ_line = (0.5 .* (Hϕ_full[:, idx_z-1] .+ Hϕ_full[:, idx_z]))[1:end-1]

    nh   = min(length(Eρ_line), grid.idx_PMLρ)
    ρ    = grid.ρd[1:nh]
    Eρ_s, Eϕ_s = Eρ_line[1:nh], Eϕ_line[1:nh]
    Hρ_s, Hϕ_s = Hρ_line[1:nh], Hϕ_line[1:nh]
    w_vec = 0.5 .* (grid.ρp[2:nh+1].^2 .- grid.ρp[1:nh].^2)

    θs = range(0, π/2, length=resolution)
    farfield_intensity = zeros(Float64, resolution)
    for (iθ, θ) in enumerate(θs)
        kρ = k0 * sin(θ); cosθ = cos(θ)
        a_plus = a_minus = b_plus = b_minus = 0.0 + 0.0im
        for i in 1:nh
            Jmp1 = besselj(m + 1, kρ * ρ[i]); Jmm1 = besselj(m - 1, kρ * ρ[i]); w = w_vec[i]
            b_plus  += ( Eϕ_s[i] + 1im*Eρ_s[i]) * Jmp1 * w
            b_minus += ( Eϕ_s[i] - 1im*Eρ_s[i]) * Jmm1 * w
            a_plus  += (-Hϕ_s[i] - 1im*Hρ_s[i]) * Jmp1 * w
            a_minus += (-Hϕ_s[i] + 1im*Hρ_s[i]) * Jmm1 * w
        end
        Eθ_ff = cosθ * ((b_minus - b_plus) + Z0*(a_minus - a_plus))
        Eϕ_ff = 1im  * ((b_plus + b_minus) + Z0*(a_plus + a_minus))
        farfield_intensity[iθ] = abs2(Eθ_ff) + abs2(Eϕ_ff)
    end
    return θs, farfield_intensity
end

function calc_flux(grid, mode, Ce;
                   flux_mon_top=nothing, flux_mon_right=nothing, flux_mon_bot=nothing,
                   report=true)
    _, Hρ, Hϕ, Hz = calc_H_fields(grid, mode, Ce)
    avg(v) = 0.5 .* (v[1:end-1] .+ v[2:end])
    mon = default_monitors(grid)

    idx_top = argmin(abs.(grid.zp .- something(flux_mon_top,   mon.top)))
    idx_ρ   = argmin(abs.(grid.ρd .- something(flux_mon_right, mon.right)))
    idx_bot = argmin(abs.(grid.zp .- something(flux_mon_bot,   mon.bot)))

    function integrate_vertical_flux(z_idx, max_r_idx)
        Eρ_line = mode.Eρ[1:max_r_idx, z_idx]
        Eϕ_line = avg(mode.Eϕ[1:max_r_idx+1, z_idx])
        Hϕ_line = 0.5 .* (Hϕ[1:max_r_idx, z_idx-1] .+ Hϕ[1:max_r_idx, z_idx])
        Hρ_line = avg(0.5 .* (Hρ[:, z_idx-1] .+ Hρ[:, z_idx])[1:max_r_idx+1])
        Sz = 0.5 .* real.(Eρ_line .* conj.(Hϕ_line) .- Eϕ_line .* conj.(Hρ_line))
        flux = 0.0
        for i in 1:max_r_idx
            flux += Sz[i] * 2π * grid.ρd[i] * grid.Δρd[i]
        end
        return flux
    end

    P_top = integrate_vertical_flux(idx_top, idx_ρ - 1)
    P_bot = integrate_vertical_flux(idx_bot, idx_ρ - 1)

    R_wall = grid.ρp[idx_ρ]; P_ρ = 0.0
    for z_idx in idx_bot:idx_top
        Ez_val = mode.Ez[idx_ρ, z_idx]
        Eϕ_val = 0.5 * (mode.Eϕ[idx_ρ, z_idx] + mode.Eϕ[idx_ρ, min(grid.Nz, z_idx+1)])
        Hz_val = 0.25 * (Hz[idx_ρ-1, z_idx] + Hz[idx_ρ, z_idx] +
                         Hz[idx_ρ-1, min(grid.Nz, z_idx+1)] + Hz[idx_ρ, min(grid.Nz, z_idx+1)])
        Hϕ_val = 0.5 * (Hϕ[idx_ρ-1, z_idx] + Hϕ[idx_ρ, z_idx])
        Sρ = 0.5 * real(Eϕ_val * conj(Hz_val) - Ez_val * conj(Hϕ_val))
        P_ρ += Sρ * grid.Δzd[z_idx] * 2π * R_wall
    end

    P_total = P_top - P_bot + P_ρ
    if report
        println("Monitor Box → R: $(round(R_wall,digits=3)) μm, Z: ",
                "$(round(grid.zp[idx_bot],digits=3)) to $(round(grid.zp[idx_top],digits=3)) μm")
        println("\n=== Power Budget ===")
        println("  Top (Up):      $(round(P_top, sigdigits=6))  ($(round(100*P_top/P_total, digits=1))%)")
        println("  Side (Out):    $(round(P_ρ,   sigdigits=6))  ($(round(100*P_ρ/P_total,   digits=1))%)")
        println("  Bottom (Down): $(round(P_bot, sigdigits=6))  ($(round(-100*P_bot/P_total,digits=1))%)")
        println("  ────────────────")
        println("  Total:         $(round(P_total, sigdigits=6))")
    end
    return P_total, P_top, P_ρ, P_bot
end

function calc_Ceff_curve(θs, farfield_intensity)
    Δθ = θs[2] - θs[1]
    Pθ = cumsum(farfield_intensity .* sin.(θs)) .* Δθ
    return Pθ ./ Pθ[end]
end

function calc_avg_emission_angle(θs, farfield_intensity)
    Δθ = θs[2] - θs[1]
    pd = farfield_intensity .* sin.(θs)
    return (sum(θs .* pd) * Δθ / (sum(pd) * Δθ)) / (π/2)
end

function calc_gaussicity(θs, farfield_intensity; NA)
    θ_w   = asin(NA)
    G     = exp.(-2.0 .* (θs ./ θ_w).^2)
    sin_θ = sin.(θs); Δθ = θs[2] - θs[1]
    overlap = abs2(sum(sqrt.(farfield_intensity) .* G .* sin_θ) * Δθ)
    norm_ff = sum(farfield_intensity .* sin_θ) * Δθ
    norm_G  = sum(G.^2 .* sin_θ) * Δθ
    return overlap / (norm_ff * norm_G)
end

default_monitors(grid) = (;
    top   = ((grid.Z - grid.PMLz_top) + grid.cbg_top) / 2,
    right = (grid.total_cbg_radius + (grid.R - grid.PMLρ)) / 2,
    bot   = (grid.PMLz_bot + grid.cbg_bot) / 2,
    farfield = ((grid.Z - grid.PMLz_top) + grid.cbg_top) / 2)

function analyze_collection(mode, grid, Ce; NA, monitors=default_monitors(grid))
    θs, I  = calc_farfield(grid, mode, Ce, monitors.farfield)
    Ceff    = calc_Ceff_curve(θs, I)
    idx_NA  = findfirst(θs .>= asin(NA))
    Ceff_NA = idx_NA === nothing ? Ceff[end] : Ceff[idx_NA]

    P_total, P_top, P_ρ, P_bot = calc_flux(grid, mode, Ce;
        flux_mon_top=monitors.top, flux_mon_right=monitors.right,
        flux_mon_bot=monitors.bot, report=false)

    η  = Ceff_NA * (P_top / P_total)
    return (; λ=mode.λ, Q=mode.Q, η, gaussicity=calc_gaussicity(θs, I; NA),
              avg_angle=calc_avg_emission_angle(θs, I),
              Ceff_NA, P_total, P_top, P_ρ, P_bot, NA, θs, farfield=I)
end

function plot_efield(grid, mode; log_scale=true, title=nothing, label=nothing,
                     x_limits=nothing, y_limits=nothing)
    xl = isnothing(x_limits) ? (0, grid.R) : x_limits
    yl = isnothing(y_limits) ? (0, grid.Z) : y_limits
    data, base = log_scale ? (log10.(mode.normE .+ 1e-20), "log₁₀|E|") :
                             (mode.normE,                  "|E|")
    default_ttl = isnothing(label) ? "$base — Electric Field" : "$label — $base"
    ttl = something(title, default_ttl)
    p = heatmap(grid.ρd, grid.zd, data',
                aspect_ratio=:equal, xlabel="ρ (μm)", ylabel="z (μm)",
                title=ttl, legend=true, xlims=xl, ylims=yl, color=:jet)
    display(p); return p
end

function plot_farfield(θs, farfield_intensity; title="Radiation Pattern", do_display=true)
    I_norm = farfield_intensity ./ maximum(farfield_intensity)
    θ_full = vcat(-reverse(θs), θs) .+ π/2
    I_full = vcat(reverse(I_norm), I_norm)
    p = plot(θ_full, I_full; proj=:polar, title, legend=false, grid=true, lw=3)
    do_display && display(p); return p
end

function plot_kspace(θs, farfield_intensity, NA; resolution=1080, do_display=true)
    kx = range(-1, 1, length=resolution); ky = range(-1, 1, length=resolution)
    img = zeros(Float64, resolution, resolution)
    for (ix, kxv) in enumerate(kx), (iy, kyv) in enumerate(ky)
        kr = hypot(kxv, kyv); kr > 1.0 && continue
        θ = asin(kr); idx = searchsortedlast(θs, θ)
        img[iy, ix] = idx == 0         ? farfield_intensity[1]   :
                      idx >= length(θs) ? farfield_intensity[end] :
                      let t = (θ - θs[idx]) / (θs[idx+1] - θs[idx])
                          (1-t)*farfield_intensity[idx] + t*farfield_intensity[idx+1]
                      end
    end
    img ./= maximum(img)
    φ = range(0, 2π, length=300)
    p = heatmap(kx, ky, img; aspect_ratio=:equal, xlabel="kx / k₀", ylabel="ky / k₀",
                title="Far-field (k-space)", color=:hot, xlims=(-1,1), ylims=(-1,1))
    plot!(p, NA.*cos.(φ), NA.*sin.(φ); color=:cyan, lw=2, linestyle=:dash, label="NA=$(NA)")
    do_display && display(p); return p
end

function plot_farfield_kspace(θs, farfield_intensity, NA; figsize=(1100, 500))
    pf = plot_farfield(θs, farfield_intensity; do_display=false)
    pk = plot_kspace(θs, farfield_intensity, NA; do_display=false)
    p  = plot(pf, pk; layout=(1,2), size=figsize)
    display(p); return p
end