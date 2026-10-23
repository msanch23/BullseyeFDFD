struct FOMTerm
    observable::Symbol
    kind::Symbol
    g::Function
    label::String
end

_lorentz_gate(o0, σ) = o -> (Δ = o - o0;      G = 1/(1 + (Δ/σ)^2); (G, -2Δ/σ^2 * G^2))
_sq_penalty(o0, w)   = o -> (Δ = (o - o0)/o0;                       (w*Δ^2, 2w*Δ/o0))
_sq_dev(o0, σ)       = o -> (Δ = (o - o0)/σ;                        (Δ^2, 2Δ/σ))

λ_target(λ0; σ=0.010, kind=:gate) =
    kind === :gate    ? FOMTerm(:λ, :gate,    _lorentz_gate(λ0, σ), "λ_target($λ0; σ=$σ)") :
    kind === :penalty ? FOMTerm(:λ, :penalty, _sq_dev(λ0, σ),       "λ_target($λ0; σ=$σ, kind=:penalty)") :
    throw(ArgumentError("λ_target: kind must be :gate or :penalty, got :$kind"))

Q_target(spec; w=1.0) =
    spec === :max ? FOMTerm(:Q, :reward,  o -> (o, 1.0),               "Q_target(:max)") :
                    FOMTerm(:Q, :penalty, _sq_penalty(float(spec), w), "Q_target($(float(spec)); w=$w)")

Fp_target(spec; w=1.0) =
    spec === :max ? FOMTerm(:Fp, :reward,  o -> (o, 1.0),               "Fp_target(:max)") :
                    FOMTerm(:Fp, :penalty, _sq_penalty(float(spec), w), "Fp_target($(float(spec)); w=$w)")

V_target(spec; w=1.0) =
    spec === :min ? FOMTerm(:V, :reward,  o -> (1/o, -1/o^2),           "V_target(:min)") :
                    FOMTerm(:V, :penalty, _sq_penalty(float(spec), w), "V_target($(float(spec)); w=$w)")

ηG_target(spec; w=1.0) =
    spec === :max ? FOMTerm(:ηG, :reward,  o -> (o, 1.0),               "ηG_target(:max)") :
                    FOMTerm(:ηG, :penalty, _sq_penalty(float(spec), w), "ηG_target($(float(spec)); w=$w)")

η_target(spec; w=1.0) =
    spec === :max ? FOMTerm(:η, :reward,  o -> (o, 1.0),               "η_target(:max)") :
                    FOMTerm(:η, :penalty, _sq_penalty(float(spec), w), "η_target($(float(spec)); w=$w)")

function evaluate_fom(terms, ovals)
    dJ_do = Dict{Symbol,Float64}()

    mult = [t for t in terms if t.kind !== :penalty]
    mv   = [t.g(ovals[t.observable]) for t in mult]
    for (i, t) in enumerate(mult)
        others = prod((first(mv[j]) for j in eachindex(mv) if j != i); init=1.0)
        dJ_do[t.observable] = get(dJ_do, t.observable, 0.0) - mv[i][2] * others
    end
    J = isempty(mv) ? 0.0 : -prod(first.(mv))

    for t in terms
        t.kind === :penalty || continue
        cost, dcost = t.g(ovals[t.observable])
        J += cost
        dJ_do[t.observable] = get(dJ_do, t.observable, 0.0) + dcost
    end
    return J, dJ_do
end

function _thickness_spec(thickness, substrates)
    map(collect(keys(thickness))) do s
        s === :tCBG && return (s, 0)
        r = match(r"^substrate(\d+)$", string(s))
        i = r === nothing ? 0 : parse(Int, r[1])
        1 <= i < length(substrates) ||
            error("thickness: :$s is not :tCBG or a finite :substrateN (N < $(length(substrates)))")
        (s, i)
    end
end

function _unpack(p, nw, tspec, substrates, tCBG)
    subs = collect(substrates); tc = tCBG
    for (k, (s, i)) in enumerate(tspec)
        s === :tCBG ? (tc = p[nw+k]) : (subs[i] = (; subs[i]..., height = p[nw+k]))
    end
    return p[1:nw], subs, tc
end

function _adjoint_setup(seed, substrates, nClad, nCBG, tCBG, slotted, λ, m, grid_kwargs; tspec=())
    nw     = length(seed) - length(tspec)
    geomof = p -> ((w, subs, tc) = _unpack(p, nw, tspec, substrates, tCBG);
                   build_geometry(w, subs, nClad, nCBG, tc; slotted))
    geom0  = geomof(seed)
    grid   = conformal_grid(geom0, λ; grid_kwargs...)
    Ce, Ch = CeChOperators(grid, λ; m=m)
    _, inv_μr = material_tensors(build_epsilon(geom0, grid), grid)
    M   = Ch * inv_μr * Ce
    εof = p -> build_epsilon(geomof(p), grid)
    return grid, M, εof, Ce, _reciprocity_weight(grid, λ; m)
end

_shifted_lu(A, k_sq) = lu(A - (k_sq + 1e-6)*I)

function _reciprocity_weight(grid, λ; m=1)
    Sρp, Sρd, Szp, Szd, ρ̃p, ρ̃d = s_vectors(grid, λ)
    aρ = ρ̃d .* grid.Δρp .* Sρd
    aϕ = ρ̃p .* grid.Δρd .* Sρp
    az = copy(aϕ); m == 0 && (az[1] = grid.ρd[1]^2 / 2)
    ℓzp = grid.Δzd .* Szp; ℓzd = grid.Δzp .* Szd
    return vcat(vec(aρ * transpose(ℓzp)), vec(-aϕ * transpose(ℓzp)), vec(az * transpose(ℓzd)))
end

function _left_eigvec(W, A, inv_ϵr, x, k_sq, grid)
    z  = W .* x ./ Vector(diag(inv_ϵr))
    ax = [2grid.N + (j-1)*grid.Nρ + 1 for j in 1:grid.Nz]
    z[ax] .= (transpose(A) * z)[ax] ./ k_sq
    y  = conj.(z) ./ norm(z)
    res = norm(A' * y .- conj(k_sq) .* y) / abs(k_sq)
    res < 1e-5 || error(@sprintf("left-eigenvector residual %.1e: W out of sync with CeChOperators", res))
    return y
end

function _adj_field(F, A, k_sq, ψ, φ, src; refine=2)
    defl!(r) = (r .-= ψ .* (dot(ψ, r)/dot(ψ, ψ)); r)
    rhs = defl!(Vector(conj.(src)))
    M   = A' - conj(k_sq)*I
    z   = F' \ ((abs(k_sq)/dot(φ, ψ)) .* ψ)
    sm(r) = (w = F' \ r; w .- z .* (dot(φ, w)/(1 + dot(φ, z))))
    v = defl!(sm(rhs))
    for _ in 1:refine
        v .+= defl!(sm(defl!(rhs .- M*v)))
    end
    return defl!(v)
end

_dk2_deps(dinv, inv_ϵr) = dinv .* .-(Vector(diag(inv_ϵr)).^2)

function _dk2_dparams(design, grad_ε, εof, grid; δ=1e-7)
    grad = zeros(ComplexF64, length(design))
    for k in eachindex(design)
        dp = copy(design); dp[k] += δ
        dm = copy(design); dm[k] -= δ
        grad[k] = sum(grad_ε .* ((_yee_eps(εof(dp), grid) .- _yee_eps(εof(dm), grid)) ./ (2δ)))
    end
    return grad
end

function _dk2_dp(design, mode, grid, A, inv_ϵr, M, εof, W; δ=1e-7)
    x = mode.raw_E
    y = _left_eigvec(W, A, inv_ϵr, x, mode.k_sq, grid)
    b = M * x
    dk2_dinv = conj(y) .* b ./ dot(y, x)
    grad_ε   = _dk2_deps(dk2_dinv, inv_ϵr)
    return _dk2_dparams(design, grad_ε, εof, grid; δ=δ), y, b
end

function _V_terms(mode, grid, ϵ, inv_ϵr; top=grid.cbg_top)
    W  = _energy_density(mode, grid, inv_ϵr)
    ir = findlast(grid.ρd .<= grid.total_cbg_radius)
    fz = [clamp((min(grid.zp[j+1], top) - max(grid.zp[j], grid.cbg_bot)) / grid.Δzp[j], 0.0, 1.0)
          for j in 1:grid.Nz]
    fz[fz .< 1e-9] .= 0.0
    ΔV  = [cell_dV(grid, i, j) for i in 1:grid.Nρ, j in 1:grid.Nz]
    box = ((1:grid.Nρ) .<= ir) .* transpose(fz)
    Wmax, istar = findmax(ifelse.(box .> 0, W, -Inf))
    Vphys = sum(W .* ΔV .* box) / Wmax
    n_cav = sqrt(maximum(real.(_cell(ϵ)[1, fz .> 0])))
    Veff  = Vphys / (mode.λ/n_cav)^3
    return (; ΔV, box, Wmax, istar, Vphys, Veff, n_cav)
end

function _dV_sources(mode, t, inv_ϵr, grid)
    N = grid.N; d = diag(inv_ϵr)
    εr = real.(reshape(1.0 ./ d[1:N],    grid.Nρ, grid.Nz))
    εϕ = real.(reshape(1.0 ./ d[N+1:2N], grid.Nρ, grid.Nz))
    εz = real.(reshape(1.0 ./ d[2N+1:3N],grid.Nρ, grid.Nz))
    κ = t.box .* t.ΔV; κ[t.istar] -= t.Vphys; κ ./= t.Wmax
    Eρ, Eϕ, Ez = mode.Eρ, mode.Eϕ, mode.Ez
    sρ = εr .* conj.(Eρ) .* κ; sϕ = εϕ .* conj.(Eϕ) .* κ; sz = εz .* conj.(Ez) .* κ
    sρ[1, :] .= (εr[1, :] .+ εϕ[1, :]) .* conj.(Eρ[1, :]) .* κ[1, :]
    sϕ[1, :] .= 0; sz[1, :] .= 0
    dρ = abs2.(Eρ) .* κ; dϕ = abs2.(Eϕ) .* κ; dz = abs2.(Ez) .* κ
    dϕ[1, :] .= abs2.(Eρ[1, :]) .* κ[1, :]; dz[1, :] .= 0
    return vcat(vec(sρ), vec(sϕ), vec(sz)), vcat(vec(dρ), vec(dϕ), vec(dz))
end

function _dV_dp(design, mode, grid, ϵ, inv_ϵr, εof, A, F, ψ, b, dk2; δ=1e-7,
                top=grid.cbg_top, it=nothing)
    x = mode.raw_E
    t = _V_terms(mode, grid, ϵ, inv_ϵr; top)
    src, dir = _dV_sources(mode, t, inv_ϵr, grid)
    v   = _adj_field(F, A, mode.k_sq, ψ, x, src)
    dAφ = _dk2_dparams(design, _dk2_deps(conj(v) .* b, inv_ϵr), εof, grid; δ=δ)
    dVp = 2 .* real.(dot(v, x) .* dk2 .- dAφ) .+
          real.(_dk2_dparams(design, dir, εof, grid; δ=δ))
    it === nothing || (dVp[it] += (_V_terms(mode, grid, ϵ, inv_ϵr; top=top+δ).Vphys -
                                   _V_terms(mode, grid, ϵ, inv_ϵr; top=top-δ).Vphys) / 2δ)
    k  = sqrt(mode.k_sq)
    dλ = -2π/real(k)^2 .* real.(dk2 ./ (2k))
    ∇V = dVp ./ (mode.λ/t.n_cav)^3 .- 3*t.Veff/mode.λ .* dλ
    return ∇V, t.Veff
end

function _farfield_terms(mode, grid, Ce; NA, m=1, ff_z=nothing, res=1080)
    N = grid.N; Z0 = 376.73; k0 = 2π/mode.λ
    raw_H = (1im/Z0) .* (Ce*mode.raw_E) ./ k0
    Hρ = reshape(raw_H[1:N], grid.Nρ, grid.Nz); Hϕ = reshape(raw_H[N+1:2N], grid.Nρ, grid.Nz)
    z0 = something(ff_z, ((grid.Z-grid.PMLz_top)+grid.cbg_top)/2); iz = argmin(abs.(grid.zp .- z0))
    avg(v) = 0.5 .* (v[1:end-1] .+ v[2:end])
    Eρl = mode.Eρ[1:end-1, iz]; Eϕl = avg(mode.Eϕ[:, iz])
    Hρl = avg(0.5 .* (Hρ[:, iz-1] .+ Hρ[:, iz])); Hϕl = (0.5 .* (Hϕ[:, iz-1] .+ Hϕ[:, iz]))[1:end-1]
    nh = min(length(Eρl), grid.idx_PMLρ); ρ = grid.ρd[1:nh]
    Eρs, Eϕs = Eρl[1:nh], Eϕl[1:nh]; Hρs, Hϕs = Hρl[1:nh], Hϕl[1:nh]
    w = 0.5 .* (grid.ρp[2:nh+1].^2 .- grid.ρp[1:nh].^2)
    θs = range(0, π/2, length=res); Δθ = θs[2]-θs[1]
    Eθ = zeros(ComplexF64, res); Eϕf = zeros(ComplexF64, res)
    for (j, θ) in enumerate(θs)
        kρ = k0*sin(θ); cθ = cos(θ); ap = am = bp = bm = 0.0+0im
        for i in 1:nh
            Jp = besselj(m+1, kρ*ρ[i]); Jm = besselj(m-1, kρ*ρ[i]); wi = w[i]
            bp += ( Eϕs[i]+1im*Eρs[i])*Jp*wi; bm += ( Eϕs[i]-1im*Eρs[i])*Jm*wi
            ap += (-Hϕs[i]-1im*Hρs[i])*Jp*wi; am += (-Hϕs[i]+1im*Hρs[i])*Jm*wi
        end
        Eθ[j] = 1im*(bp+bm) - Z0*cθ*(am-ap); Eϕf[j] = cθ*(bm-bp) + 1im*Z0*(ap+am)
    end
    I = abs2.(Eθ) .+ abs2.(Eϕf); sinθ = sin.(θs)
    θw = asin(NA); G = exp.(-2.0 .* (θs./θw).^2); sqrtI = sqrt.(I)
    S = sum(sqrtI.*G.*sinθ)*Δθ; P = sum(I.*sinθ)*Δθ; normG = sum(G.^2 .*sinθ)*Δθ
    ηG = S^2/(P*normG)
    Pθ = cumsum(I.*sinθ).*Δθ; idxNA = findfirst(θs .>= θw)
    Ceff = idxNA === nothing ? 1.0 : Pθ[idxNA]/Pθ[end]
    return (; θs, Δθ, sinθ, G, sqrtI, Eθ, Eϕf, S, P, normG, ηG, Ceff, idxNA, ρ, w, nh, iz, k0)
end

function _ff_backward(ft, grid, Ce, dI; m=1)
    N = grid.N; Z0 = 376.73; k0 = ft.k0; nh = ft.nh; iz = ft.iz; ρ = ft.ρ; w = ft.w
    dEθ = dI .* conj.(ft.Eθ); dEϕ = dI .* conj.(ft.Eϕf)
    dEρs = zeros(ComplexF64,nh); dEϕs = zeros(ComplexF64,nh)
    dHρs = zeros(ComplexF64,nh); dHϕs = zeros(ComplexF64,nh)
    for (j, θ) in enumerate(ft.θs)
        kρ = k0*sin(θ); cθ = cos(θ); gθ = dEθ[j]; gϕ = dEϕ[j]
        dbm = gθ*1im+gϕ*cθ; dbp = gθ*1im-gϕ*cθ; dam = -gθ*cθ*Z0+gϕ*1im*Z0; dap = gθ*cθ*Z0+gϕ*1im*Z0
        for i in 1:nh
            Jp = besselj(m+1, kρ*ρ[i]); Jm = besselj(m-1, kρ*ρ[i]); wi = w[i]
            dEϕs[i] += dbp*Jp*wi + dbm*Jm*wi
            dEρs[i] += dbp*1im*Jp*wi + dbm*(-1im)*Jm*wi
            dHϕs[i] += dap*(-1.0)*Jp*wi + dam*(-1.0)*Jm*wi
            dHρs[i] += dap*(-1im)*Jp*wi + dam*(1im)*Jm*wi
        end
    end
    dEρ = zeros(ComplexF64,grid.Nρ,grid.Nz); dEϕ = zeros(ComplexF64,grid.Nρ,grid.Nz)
    dHρ = zeros(ComplexF64,grid.Nρ,grid.Nz); dHϕ = zeros(ComplexF64,grid.Nρ,grid.Nz)
    for i in 1:nh
        dEρ[i,iz] += dEρs[i]; dEϕ[i,iz] += 0.5*dEϕs[i]; dEϕ[i+1,iz] += 0.5*dEϕs[i]
        dHρ[i,iz-1] += 0.25*dHρs[i]; dHρ[i,iz] += 0.25*dHρs[i]
        dHρ[i+1,iz-1] += 0.25*dHρs[i]; dHρ[i+1,iz] += 0.25*dHρs[i]
        dHϕ[i,iz-1] += 0.5*dHϕs[i]; dHϕ[i,iz] += 0.5*dHϕs[i]
    end
    dH = zeros(ComplexF64,3N); dH[1:N] .= vec(dHρ); dH[N+1:2N] .= vec(dHϕ)
    dφE = zeros(ComplexF64,3N); dφE[1:N] .= vec(dEρ); dφE[N+1:2N] .= vec(dEϕ)
    return dφE .+ (1im/(Z0*k0)) .* (transpose(Ce)*dH)
end

function _flux_terms(mode, grid, Ce)
    N = grid.N; Z0 = 376.73; k0 = 2π/mode.λ
    raw_H = (1im/Z0) .* (Ce*mode.raw_E) ./ k0
    Hρ = reshape(raw_H[1:N],grid.Nρ,grid.Nz); Hϕ = reshape(raw_H[N+1:2N],grid.Nρ,grid.Nz)
    Hz = reshape(raw_H[2N+1:3N],grid.Nρ,grid.Nz)
    Eρ, Eϕ, Ez = mode.Eρ, mode.Eϕ, mode.Ez
    mon = default_monitors(grid)
    itop = argmin(abs.(grid.zp .- mon.top)); iρ = argmin(abs.(grid.ρd .- mon.right))
    ibot = argmin(abs.(grid.zp .- mon.bot)); max_r = iρ-1; avg(v) = 0.5 .* (v[1:end-1] .+ v[2:end])
    function vflux(z)
        Eρl = Eρ[1:max_r,z]; Eϕl = avg(Eϕ[1:max_r+1,z])
        Hϕl = 0.5 .* (Hϕ[1:max_r,z-1] .+ Hϕ[1:max_r,z]); Hρl = avg(0.5 .* (Hρ[:,z-1] .+ Hρ[:,z])[1:max_r+1])
        Sz = 0.5 .* real.(Eρl.*conj.(Hϕl) .- Eϕl.*conj.(Hρl))
        s = 0.0; for i in 1:max_r; s += Sz[i]*π*(grid.ρp[i+1]^2-grid.ρp[i]^2); end; s
    end
    P_top = vflux(itop); P_bot = vflux(ibot); R_wall = grid.ρp[iρ]; P_ρ = 0.0
    for z in ibot:itop-1
        zn = min(grid.Nz,z+1); Ez_v = Ez[iρ,z]; Eϕ_v = 0.5*(Eϕ[iρ,z]+Eϕ[iρ,zn])
        Hz_v = 0.25*(Hz[iρ-1,z]+Hz[iρ,z]+Hz[iρ-1,zn]+Hz[iρ,zn]); Hϕ_v = 0.5*(Hϕ[iρ-1,z]+Hϕ[iρ,z])
        P_ρ += 0.5*real(Eϕ_v*conj(Hz_v)-Ez_v*conj(Hϕ_v)) * grid.Δzp[z]*2π*R_wall
    end
    P_total = P_top-P_bot+P_ρ
    return (; P_top, P_bot, P_ρ, P_total, top_frac=P_top/P_total,
              itop, ibot, iρ, max_r, R_wall, k0, Eρ, Eϕ, Ez, Hρ, Hϕ, Hz)
end

function _dtopfrac_de(fx, grid, Ce)
    N = grid.N; Z0 = 376.73; k0 = fx.k0; tf = fx.top_frac; Pt = fx.P_total
    max_r = fx.max_r; iρ = fx.iρ; avg(v) = 0.5 .* (v[1:end-1] .+ v[2:end])
    dEρ = zeros(ComplexF64,grid.Nρ,grid.Nz); dEϕ = zeros(ComplexF64,grid.Nρ,grid.Nz); dEz = zeros(ComplexF64,grid.Nρ,grid.Nz)
    dHρ = zeros(ComplexF64,grid.Nρ,grid.Nz); dHϕ = zeros(ComplexF64,grid.Nρ,grid.Nz); dHz = zeros(ComplexF64,grid.Nρ,grid.Nz)
    for (z, dfdP) in ((fx.itop, (1-tf)/Pt), (fx.ibot, tf/Pt))
        Eρl = fx.Eρ[1:max_r,z]; Eϕl = avg(fx.Eϕ[1:max_r+1,z])
        Hϕl = 0.5 .* (fx.Hϕ[1:max_r,z-1] .+ fx.Hϕ[1:max_r,z]); Hρl = avg(0.5 .* (fx.Hρ[:,z-1] .+ fx.Hρ[:,z])[1:max_r+1])
        for i in 1:max_r
            coeff = dfdP*π*(grid.ρp[i+1]^2-grid.ρp[i]^2)
            dEρv = 0.25*conj(Hϕl[i])*coeff; dEϕv = -0.25*conj(Hρl[i])*coeff
            dHϕv = 0.25*conj(Eρl[i])*coeff; dHρv = -0.25*conj(Eϕl[i])*coeff
            dEρ[i,z] += dEρv; dEϕ[i,z] += 0.5*dEϕv; dEϕ[i+1,z] += 0.5*dEϕv
            dHϕ[i,z-1] += 0.5*dHϕv; dHϕ[i,z] += 0.5*dHϕv
            dHρ[i,z-1] += 0.25*dHρv; dHρ[i,z] += 0.25*dHρv; dHρ[i+1,z-1] += 0.25*dHρv; dHρ[i+1,z] += 0.25*dHρv
        end
    end
    for z in fx.ibot:fx.itop-1
        zn = min(grid.Nz,z+1); coeff = (-tf/Pt)*grid.Δzp[z]*2π*fx.R_wall
        Eϕ_v = 0.5*(fx.Eϕ[iρ,z]+fx.Eϕ[iρ,zn]); Ez_v = fx.Ez[iρ,z]
        Hz_v = 0.25*(fx.Hz[iρ-1,z]+fx.Hz[iρ,z]+fx.Hz[iρ-1,zn]+fx.Hz[iρ,zn]); Hϕ_v = 0.5*(fx.Hϕ[iρ-1,z]+fx.Hϕ[iρ,z])
        dEϕv = 0.25*conj(Hz_v)*coeff; dEzv = -0.25*conj(Hϕ_v)*coeff
        dHzv = 0.25*conj(Eϕ_v)*coeff; dHϕv = -0.25*conj(Ez_v)*coeff
        dEϕ[iρ,z] += 0.5*dEϕv; dEϕ[iρ,zn] += 0.5*dEϕv; dEz[iρ,z] += dEzv
        dHz[iρ-1,z] += 0.25*dHzv; dHz[iρ,z] += 0.25*dHzv; dHz[iρ-1,zn] += 0.25*dHzv; dHz[iρ,zn] += 0.25*dHzv
        dHϕ[iρ-1,z] += 0.5*dHϕv; dHϕ[iρ,z] += 0.5*dHϕv
    end
    dH = zeros(ComplexF64,3N); dH[1:N] .= vec(dHρ); dH[N+1:2N] .= vec(dHϕ); dH[2N+1:3N] .= vec(dHz)
    dφE = zeros(ComplexF64,3N); dφE[1:N] .= vec(dEρ); dφE[N+1:2N] .= vec(dEϕ); dφE[2N+1:3N] .= vec(dEz)
    return dφE .+ (1im/(Z0*k0)) .* (transpose(Ce)*dH)
end

function _field_grad(design, x, k_sq, grid, inv_ϵr, εof, A, F, ψ, b, dk2, src, gλ; δ=1e-7)
    v   = _adj_field(F, A, k_sq, ψ, x, src)
    dAφ = _dk2_dparams(design, _dk2_deps(conj(v) .* b, inv_ϵr), εof, grid; δ=δ)
    eig = 2 .* real.(dot(v, x) .* dk2 .- dAφ)
    k   = sqrt(k_sq); dλ = -2π/real(k)^2 .* real.(dk2 ./ (2k))
    return eig .+ gλ .* dλ
end

function _dηG_dp(design, mode, grid, inv_ϵr, εof, A, F, ψ, b, dk2, Ce; NA, δ=1e-7, δλ=1e-7)
    ft  = _farfield_terms(mode, grid, Ce; NA)
    dI  = (ft.sinθ .* ft.Δθ ./ ft.P) .* (ft.S .* ft.G ./ ((ft.sqrtI .+ 1e-30) .* ft.normG) .- ft.ηG)
    gλ  = (_farfield_terms((; mode..., λ=mode.λ+δλ), grid, Ce; NA).ηG -
           _farfield_terms((; mode..., λ=mode.λ-δλ), grid, Ce; NA).ηG) / (2δλ)
    ∇ηG = _field_grad(design, mode.raw_E, mode.k_sq, grid, inv_ϵr, εof, A, F, ψ, b, dk2,
                      _ff_backward(ft, grid, Ce, dI), gλ; δ=δ)
    return ∇ηG, ft.ηG
end

function _dη_dp(design, mode, grid, inv_ϵr, εof, A, F, ψ, b, dk2, Ce; NA, δ=1e-7, δλ=1e-7)
    ft = _farfield_terms(mode, grid, Ce; NA); fx = _flux_terms(mode, grid, Ce)
    Ceff = ft.Ceff; tf = fx.top_frac
    dCeff_dI = (ft.Δθ ./ ft.P) .* (((1:length(ft.θs)) .<= ft.idxNA) .- Ceff) .* ft.sinθ
    src = tf .* _ff_backward(ft, grid, Ce, dCeff_dI) .+ Ceff .* _dtopfrac_de(fx, grid, Ce)
    dCeff_dλ = (_farfield_terms((; mode..., λ=mode.λ+δλ), grid, Ce; NA).Ceff -
                _farfield_terms((; mode..., λ=mode.λ-δλ), grid, Ce; NA).Ceff) / (2δλ)
    ∇η = _field_grad(design, mode.raw_E, mode.k_sq, grid, inv_ϵr, εof, A, F, ψ, b, dk2,
                     src, tf*dCeff_dλ; δ=δ)
    return ∇η, Ceff*tf
end

function _fom_and_grad(k_sq, dk2, terms; fields=Dict{Symbol,Tuple{Float64,Vector{Float64}}}())
    k = sqrt(k_sq); kr = real(k); ki = imag(k); sgn = sign(ki)
    λ = 2π/kr; Q = kr/(2abs(ki))

    ovals = Dict{Symbol,Float64}(:λ => λ, :Q => Q)
    for (s, (val, _)) in fields; ovals[s] = val; end
    haskey(ovals, :V) && (ovals[:Fp] = (3/(4π^2)) * Q / ovals[:V])
    J, dJ_do = evaluate_fom(terms, ovals)

    if haskey(dJ_do, :Fp)
        dJ_do[:Q] = get(dJ_do, :Q, 0.0) + dJ_do[:Fp] * (3/(4π^2)) / ovals[:V]
        dJ_do[:V] = get(dJ_do, :V, 0.0) - dJ_do[:Fp] * (3/(4π^2)) * Q / ovals[:V]^2
    end

    ∇ = zeros(Float64, length(dk2))
    for i in eachindex(dk2)
        dk = dk2[i]/(2k); dkr = real(dk); dki = imag(dk)
        dλ = -2π/kr^2 * dkr
        dQ = dkr/(2abs(ki)) - Q/abs(ki)*sgn*dki
        ∇[i] = get(dJ_do, :λ, 0.0)*dλ + get(dJ_do, :Q, 0.0)*dQ
    end
    for (s, (_, ∇field)) in fields
        c = get(dJ_do, s, 0.0)
        c == 0.0 || (∇ .+= c .* ∇field)
    end
    return J, ∇, ovals
end

_fmt_obs(s, v) = s === :λ ? @sprintf("%.6f", v) :
                 s === :Q ? @sprintf("%.2f", v) :
                            @sprintf("%.6g", v)

_git(args...) = try readchomp(pipeline(`git -C $(pkgdir(@__MODULE__)) $args`; stderr=devnull)) catch; "" end

# Uncommitted src/ + Project.toml changes, so a "-dirty" run can be reproduced (commit + patch).
# Written once per run folder, plus code_round<k>.patch if a later round ran different code.
function _save_patch(rundir, k)
    d = _git("diff", "HEAD", "--", "src", "Project.toml") * "\n"
    p = joinpath(rundir, "code.patch")
    isfile(p) && read(p, String) == d && return
    write(isfile(p) ? joinpath(rundir, "code_round$k.patch") : p, d)
end

# Log layout: results/<yymmdd-HHMMSS>[_label]/round<k>.csv, one folder per run. A round seeded from a
# previous result (adj_opt(prev, …)) writes the next round<k>.csv into that result's folder. With
# save_data = a .jl path (@__FILE__) the folder sits next to that script and keeps a copy of it.
# An explicit .csv path is used as given (no run folder, no index).
function _open_log(save_data, label, prev, seed, cols, log_obs, meta)
    save_data === false && return nothing, nothing, 1
    pf = prev === nothing || !hasproperty(prev, :logfile) ? nothing : prev.logfile
    k  = pf === nothing ? 1 : (hasproperty(prev, :round) ? prev.round : 1) + 1
    if save_data isa AbstractString && endswith(save_data, ".csv")
        path = save_data; mkpath(dirname(path))
    else
        script = save_data isa AbstractString && endswith(save_data, ".jl") ? abspath(save_data) : nothing
        base   = script !== nothing ? dirname(script) : save_data isa AbstractString ? save_data :
                 (isempty(PROGRAM_FILE) ? pwd() : dirname(abspath(PROGRAM_FILE)))
        rundir = pf !== nothing && occursin(r"^round\d+\.csv$", basename(pf)) ? dirname(pf) :
                 joinpath(base, "results", Dates.format(now(), "yymmdd-HHMMSS") * (isempty(label) ? "" : "_$label"))
        mkpath(rundir)
        path = joinpath(rundir, "round$k.csv")
        isfile(path) && (path = joinpath(rundir, "round$k-" * Dates.format(now(), "HHMMSS") * ".csv"))
        script === nothing || isfile(joinpath(rundir, basename(script))) || cp(script, joinpath(rundir, basename(script)))
        endswith(_git("describe", "--always", "--dirty", "--tags"), "-dirty") && _save_patch(rundir, k)
    end
    io = open(path, "w")
    gh = _git("describe", "--always", "--dirty", "--tags")
    println(io, "# git=$(isempty(gh) ? "unknown" : gh)  version=$(pkgversion(@__MODULE__))")
    println(io, "# run: $(basename(dirname(path)))  round $k" *
                (pf === nothing ? "" : "  (seeded from $(basename(pf))" *
                                       (hasproperty(prev, :best_eval) ? " eval $(prev.best_eval))" : ")")))
    for line in meta; println(io, "# ", line); end
    println(io, "# seed=", seed)
    println(io, "eval,elapsed_s,J,", join(string.(log_obs), ","), ",", join(cols, ","))
    flush(io)
    return io, path, k
end

const _INDEX_OBS = (:λ, :Q, :V, :Fp, :η, :ηG, :overlap_seed)

# One line per finished round in results/index.csv (only for logs inside a run folder).
function _index_row(path, k, started, optimizer, thickness, grid_kwargs, np, nevals, b, stop)
    occursin(r"^round\d+(-\d+)?\.csv$", basename(path)) || return
    idx = joinpath(dirname(dirname(path)), "index.csv")
    new = !isfile(idx)
    open(idx, "a") do io
        new && println(io, "run,round,started,git,optimizer,thickness,grid,params,evals,returned_eval,J,",
                           join(string.(_INDEX_OBS), ","), ",stop")
        grid = isempty(grid_kwargs) ? "default" : join(("$s=$v" for (s, v) in pairs(grid_kwargs)), " ")
        thk  = isempty(thickness) ? "fixed" : join(string.(keys(thickness)), "+")
        obs  = join((haskey(b.ob, s) ? _fmt_obs(s, b.ob[s]) : "" for s in _INDEX_OBS), ",")
        println(io, join((basename(dirname(path)), k, Dates.format(started, "yyyy-mm-dd HH:MM"),
                          _git("describe", "--always", "--dirty", "--tags"), optimizer, thk, grid, np, nevals,
                          b.n, @sprintf("%.6f", b.J)), ","), ",", obs, ",\"", stop, "\"")
    end
end

function _log_row(io, log_obs, n_eval, elapsed, J, ovals, x)
    io === nothing && return
    obs  = join((haskey(ovals, s) ? _fmt_obs(s, ovals[s]) : "NaN" for s in log_obs), ",")
    pstr = join((@sprintf("%.6f", xi) for xi in x), ",")
    @printf(io, "%d,%.1f,%.6f,%s,%s\n", n_eval, elapsed, J, obs, pstr)
    flush(io)
end

_adj_options(iterations) = Optim.Options(; iterations, outer_iterations=10, x_abstol=1e-3,
                                           f_reltol=1e-4, outer_f_reltol=1e-4, show_trace=true)

const _FTOL     = 1e-4   # L-BFGS-B stops once each of _PATIENCE iterations in a row lowers J by ≤ _FTOL·|J|
const _PATIENCE = 3      # (relative to J itself: Fortran's own factr test floors |J| at 1, i.e. absolute for |J| < 1)
const _PGTOL    = 1e-5   # … or once the largest projected-gradient entry is below this

# L-BFGS-B (Byrd, Lu, Nocedal & Zhu): bounds enforced by projection, so a parameter can sit on its
# bound while the others keep moving; no barrier term enters J. The reverse-communication loop is
# driven here so the stopping test can be relative.
function _lbfgsb(f, g!, x0, lower, upper; maxiter, m=10)
    n = length(x0)
    o = L_BFGS_B(n, m)
    o.l .= something(lower, fill(-Inf, n)); o.u .= something(upper, fill(Inf, n))
    o.nbd .= [isinf(l) ? (isinf(u) ? 0 : 3) : (isinf(u) ? 1 : 2) for (l, u) in zip(o.l, o.u)]
    x = copy(x0); fx = Inf; Jprev = Inf; slow = 0; it = 0
    done(msg) = (; minimum = fx, minimizer = x, iterations = it, message = msg)
    o.task[1:5] = b"START"
    while true
        setulb(n, m, x, o.l, o.u, o.nbd, fx, o.g, 0.0, _PGTOL, o.wa, o.iwa, o.task, -1,
               o.csave, o.lsave, o.isave, o.dsave)
        if o.task[1:2] == b"FG"
            fx = f(x); g!(o.g, x)
        elseif o.task[1:5] == b"NEW_X"
            it += 1
            slow  = Jprev - fx ≤ _FTOL * abs(fx) ? slow + 1 : 0
            Jprev = fx
            it ≥ maxiter        && return done("STOP: iteration limit ($maxiter)")
            slow ≥ _PATIENCE    && return done("CONVERGENCE: ΔJ ≤ $(_FTOL)·|J| for $(_PATIENCE) iterations")
        else
            return done(rstrip(String(copy(o.task))))
        end
    end
end

_overlap(a, b, w) = abs(sum(w .* a .* conj.(b))) / sqrt(sum(w .* abs2.(a)) * sum(w .* abs2.(b)))

_valid_modes(raw, Qfloor, grid, inv_ϵr) = filter(m -> m.Q ≥ Qfloor && _pml_fraction(m, grid, inv_ϵr) ≤ 0.5, raw)

function _nearest_mode(raw, ref, Qfloor, grid, inv_ϵr)
    cand = _valid_modes(raw, Qfloor, grid, inv_ϵr)
    isempty(cand) && return nothing
    argmin(m -> ((m.λ - ref.λ)/0.001)^2 + ((m.Q - ref.Q)/(0.05*ref.Q))^2, cand)
end

function _follow_mode(raw, xref, w, track, Qfloor, grid, inv_ϵr)
    cand = _valid_modes(raw, Qfloor, grid, inv_ϵr)
    isempty(cand) && return nothing
    o, i = findmax(m -> _overlap(m.raw_E, xref, w), cand)
    return o ≥ track ? cand[i] : nothing
end

# Next remesh round: seed from a previous adj_opt result, which also carries the mode to follow.
function adj_opt(prev::NamedTuple, FOMs; kwargs...)
    (haskey(kwargs, :substrates) || haskey(kwargs, :tCBG)) &&
        throw(ArgumentError("adj_opt(prev, …) takes substrates and tCBG from the previous result"))
    return adj_opt(prev.design, FOMs; substrates = prev.substrates, tCBG = prev.tCBG, ref = prev, kwargs...)
end

function adj_opt(seed_design, FOMs;
                 substrates, nClad, nCBG, tCBG, slotted=false, thickness=(;),
                 λ=0.780, NA=nothing, lower=nothing, upper=nothing,
                 Nmodes=10, m=1, grid_kwargs=(;),
                 σ_sel=0.020, Qfloor=30.0, track=0.5, ref=nothing, save_data=false, label="",
                 optimizer=:lbfgsb, iterations=nothing,
                 linesearch=BackTracking(), mu0=NaN, options=nothing)   # the last three: :fminbox only

    m == 1 || error("adj_opt supports m = 1 only in this version (got m = $m)")
    optimizer in (:fminbox, :lbfgsb) ||
        throw(ArgumentError("optimizer must be :fminbox or :lbfgsb, got :$optimizer"))
    iterations = something(iterations, optimizer === :lbfgsb ? 200 : 50)   # L-BFGS-B: a cap; Fminbox: per inner run
    options    = something(options, _adj_options(iterations))
    tspec = _thickness_spec(thickness, substrates)
    nw    = length(seed_design)
    p0    = [seed_design; Float64[s === :tCBG ? tCBG : substrates[i].height for (s, i) in tspec]]
    it    = (k = findfirst(s -> s[1] === :tCBG, tspec)) === nothing ? nothing : nw + k
    if !isempty(tspec)
        lower = [something(lower, fill(-Inf, nw)); Float64[first(b) for b in values(thickness)]]
        upper = [something(upper, fill( Inf, nw)); Float64[last(b)  for b in values(thickness)]]
    end

    grid, M, εof, Ce, W = _adjoint_setup(p0, substrates, nClad, nCBG, tCBG,
                                         slotted, λ, m, grid_kwargs; tspec)
    if !isempty(tspec)
        _, subs_hi, tc_hi = _unpack(upper, nw, tspec, substrates, tCBG)
        grid.cbg_bot - sum(s.height for s in subs_hi[1:end-1]; init=0.0) > grid.PMLz_bot ||
            error("thickness upper bounds push the substrate stack into the bottom PML of the seed grid")
        grid.cbg_bot + tc_hi < default_monitors(grid).top ||
            error("tCBG upper bound reaches the top flux/far-field plane of the seed grid")
    end

    log_obs  = [unique([t.observable for t in FOMs]); :pml_frac;
                track === nothing ? Symbol[] : [:overlap, :overlap_seed]]
    refλQ    = ref === nothing ? nothing : hasproperty(ref, :mode) ? ref.mode : ref
    wv       = abs.(W)
    xseed    = Ref{Union{Nothing,Vector{ComplexF64}}}(nothing)   # mode at the seed
    xref     = Ref{Union{Nothing,Vector{ComplexF64}}}(nothing)   # mode at the best point so far
    best     = Ref{Any}(nothing)                                 # lowest-J evaluation: (; n, x, J, λ, Q)
    needs_V  = any(t.observable in (:V, :Fp) for t in FOMs)
    needs_ηG = any(t.observable === :ηG for t in FOMs)
    needs_η  = any(t.observable === :η  for t in FOMs)
    needs_F  = needs_V || needs_ηG || needs_η
    (needs_ηG || needs_η) && NA === nothing &&
        error("η/ηG figures of merit require an NA (pass NA=… to adj_opt)")
    meta = ["cbg: nClad=$nClad  nCBG=$nCBG  tCBG=$tCBG  slotted=$slotted",
            "substrates: " * join(["(n=$(s.n), h=$(s.height))" for s in substrates], "  "),
            "solve: grid_kwargs=$grid_kwargs  Nmodes=$Nmodes  m=$m  λ=$λ",
            "FOM: " * join([t.label for t in FOMs], " · ") * (NA === nothing ? "" : "  (NA=$NA)"),
            "select: Qfloor=$Qfloor  PML fraction ≤ 0.5  σ_sel=$σ_sel"]
    push!(meta, lower === nothing ? "bounds: none" : "bounds: p1–p$nw lower=$(lower[1:nw])  upper=$(upper[1:nw])")
    push!(meta, optimizer === :lbfgsb ?
                "optimizer: L-BFGS-B  m=10  maxiter=$iterations  stop: ΔJ ≤ $(_FTOL)·|J| for $(_PATIENCE) iterations, pgtol=$_PGTOL" :
                "optimizer: $(lower === nothing ? "LBFGS" : "Fminbox(LBFGS)")  iterations=$iterations" *
                "  mu0=$(isnan(mu0) ? "auto" : mu0)")
    isempty(tspec) || push!(meta, "thickness: " * join(("$s ∈ $b" for (s, b) in pairs(thickness)), "  "))
    track === nothing || push!(meta, "track: highest overlap with the best point's mode, ≥ $track" *
                                     (refλQ === nothing ? "" : "; seed mode nearest λ=$(refλQ.λ), Q=$(refλQ.Q)"))
    nsub = length(substrates) - 1                                # finite layers; the last one is the wafer
    cols = [["p$i" for i in 1:nw]; "tCBG"; ["substrate$i" for i in 1:nsub]]
    row(x) = ((w, subs, tc) = _unpack(x, nw, tspec, substrates, tCBG); [w; tc; [s.height for s in subs[1:nsub]]])
    started = now()
    logio, logpath, rnd = _open_log(save_data, label, ref, row(p0), cols, log_obs, meta)

    cache_x = Ref(Float64[]); cache = Ref{Any}(nothing); n_eval = Ref(0)
    t_start = time()

    function evaluate(x)
        if cache_x[] != x
            cache_x[] = copy(x); n_eval[] += 1
            ϵ      = εof(x)
            inv_ϵr = material_tensors(ϵ, grid)[1]
            A      = inv_ϵr * M
            raw = _solve_eig(A, grid, λ, Nmodes, 1e-6, 3000)
            rm  = if track === nothing || (xref[] === nothing && refλQ === nothing)
                      _select_mode(raw, λ, σ_sel; Qfloor, grid, inv_ϵr)
                  elseif xref[] === nothing                  # seed of a later round: the previous round's mode
                      _nearest_mode(raw, refλQ, Qfloor, grid, inv_ϵr)
                  else
                      _follow_mode(raw, xref[], wv, track, Qfloor, grid, inv_ϵr)
                  end
            if rm === nothing
                @printf("[%3d] no valid mode (none with Q≥%.0f, PML fraction ≤0.5%s)  ‖Δp‖=%.3f\n",
                        n_eval[], Qfloor, xref[] === nothing ? "" : ", overlap ≥ $track", norm(x .- p0))
                for (i, md) in enumerate(raw)
                    @printf("        raw %d: λ=%.1f nm  Q=%.1f  PML %.2f%s\n", i, md.λ*1000, md.Q,
                            _pml_fraction(md, grid, inv_ϵr),
                            xref[] === nothing ? "" : @sprintf("  overlap %.2f", _overlap(md.raw_E, xref[], wv)))
                end
                n_eval[] == 1 && error("adj_opt: no valid mode at the seed design (raw modes listed above)")
                cache[] = (J = 1e3, ∇ = zeros(length(x)))
                _log_row(logio, log_obs, n_eval[], time()-t_start, 1e3, Dict{Symbol,Float64}(), row(x))
            else
                dk2, y, b = _dk2_dp(x, rm, grid, A, inv_ϵr, M, εof, W; δ=1e-7)
                F = needs_F ? _shifted_lu(A, rm.k_sq) : nothing
                fields = Dict{Symbol,Tuple{Float64,Vector{Float64}}}()
                if needs_V
                    top = grid.cbg_bot + _unpack(x, nw, tspec, substrates, tCBG)[3]
                    ∇V, Veff = _dV_dp(x, rm, grid, ϵ, inv_ϵr, εof, A, F, y, b, dk2; δ=1e-7, top, it)
                    fields[:V] = (Veff, ∇V)
                end
                if needs_ηG
                    ∇ηG, ηGv = _dηG_dp(x, rm, grid, inv_ϵr, εof, A, F, y, b, dk2, Ce; NA)
                    fields[:ηG] = (ηGv, ∇ηG)
                end
                if needs_η
                    ∇η, ηv = _dη_dp(x, rm, grid, inv_ϵr, εof, A, F, y, b, dk2, Ce; NA)
                    fields[:η] = (ηv, ∇η)
                end
                J, ∇, ob = _fom_and_grad(rm.k_sq, dk2, FOMs; fields)
                ob[:pml_frac] = _pml_fraction(rm, grid, inv_ϵr)
                if track !== nothing
                    ob[:overlap]      = xref[]  === nothing ? 1.0 : _overlap(rm.raw_E, xref[], wv)
                    ob[:overlap_seed] = xseed[] === nothing ? 1.0 : _overlap(rm.raw_E, xseed[], wv)
                    xseed[] === nothing && (xseed[] = copy(rm.raw_E))
                end
                if best[] === nothing || J < best[].J
                    best[] = (; n = n_eval[], x = copy(x), J, λ = ob[:λ], Q = ob[:Q], ob = copy(ob))
                    track === nothing || (xref[] = copy(rm.raw_E))
                end
                cache[]  = (; J, ∇)
                @printf("[%3d] J=%.4f  λ=%.2f nm  Q=%.1f%s%s\n", n_eval[], J, ob[:λ]*1000,
                        ob[:Q], haskey(ob, :V) ? @sprintf("  V=%.4f", ob[:V]) : "",
                        track === nothing ? "" : @sprintf("  ovl=%.2f seed=%.2f", ob[:overlap], ob[:overlap_seed]))
                _log_row(logio, log_obs, n_eval[], time()-t_start, J, ob, row(x))
            end
        end
        return cache[]
    end
    f(x)     = evaluate(x).J
    g!(G, x) = (G .= evaluate(x).∇)

    result = try
        optimizer === :lbfgsb ? _lbfgsb(f, g!, p0, lower, upper; maxiter=iterations) :
        lower === nothing     ? optimize(f, g!, p0, LBFGS(; linesearch), options) :
                                optimize(f, g!, lower, upper, p0, Fminbox(LBFGS(; linesearch); mu0), options)
    finally
        logio !== nothing && close(logio)
    end

    # Return the lowest-J evaluation, not the optimizer's last iterate.
    b    = best[]
    stop = optimizer === :lbfgsb ? "L-BFGS-B: $(result.message) after $(result.iterations) iterations" :
           "$(lower === nothing ? "LBFGS" : "Fminbox"): $(Optim.converged(result) ? "converged" : "not converged")" *
           " after $(Optim.iterations(result)) iterations"
    println("adj_opt: ", stop)
    @printf("adj_opt: returning eval %d (lowest J = %.4f)\n", b.n, b.J)
    if logpath !== nothing
        open(logpath, "a") do io
            println(io, "# stop: ", stop)
            println(io, "# returned: eval $(b.n) (lowest J)")
        end
        _index_row(logpath, rnd, started, optimizer, thickness, grid_kwargs, length(p0), n_eval[], b, stop)
    end
    w, subs, tc = _unpack(b.x, nw, tspec, substrates, tCBG)
    return (; design = w, tCBG = tc, substrates = subs, fom = b.J, mode = (; λ = b.λ, Q = b.Q),
              round = rnd, best_eval = b.n, grid, result, logfile = logpath)
end

