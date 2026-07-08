_eps(n) = ComplexF64(n)^2

function build_epsilon(geometry, grid)
    ε = Matrix{ComplexF64}(undef, grid.Nρ, grid.Nz)

    N        = length(geometry.ring.width)
    has_slot = haskey(geometry, :slot)
    widths   = has_slot ? Float64[geometry.slot.width, geometry.eye.width]        : Float64[geometry.eye.width]
    eps_seg  = has_slot ? ComplexF64[_eps(geometry.slot.n), _eps(geometry.eye.n)] : ComplexF64[_eps(geometry.eye.n)]

    for i in 1:N
        push!(widths, geometry.trench.width[i]); push!(eps_seg, _eps(geometry.trench.n))
        push!(widths, geometry.ring.width[i]);   push!(eps_seg, _eps(geometry.ring.n))
    end
    push!(widths, geometry.buffer.width); push!(eps_seg, _eps(geometry.buffer.n))
    push!(eps_seg, _eps(geometry.eye.n))
    ρ_interfaces = cumsum(widths)

    function radial_blend(ρ_lo, ρ_hi)
        A = ρ_hi^2 - ρ_lo^2
        A <= 0 && return eps_seg[min(searchsortedfirst(ρ_interfaces, ρ_lo), length(eps_seg))]
        ϵ = zero(ComplexF64)
        seg = searchsortedfirst(ρ_interfaces, ρ_lo)
        while ρ_lo < ρ_hi && seg <= length(eps_seg)
            edge = seg <= length(ρ_interfaces) ? min(ρ_interfaces[seg], ρ_hi) : ρ_hi
            if edge > ρ_lo
                ϵ += ((edge^2 - ρ_lo^2) / A) * eps_seg[seg]
            end
            ρ_lo = edge
            seg += 1
        end
        return ϵ
    end

    ε_band = [radial_blend(grid.ρp[i], grid.ρp[i+1]) for i in 1:grid.Nρ]

    cbg_bot = grid.cbg_bot
    cbg_top = cbg_bot + geometry.tCBG
    zlayers = Tuple{Float64,Float64,ComplexF64}[]
    z_top   = cbg_bot
    for sub in Base.front(geometry.stack)
        z_bot = z_top - sub.height
        push!(zlayers, (z_bot, z_top, _eps(sub.n)))
        z_top = z_bot
    end
    push!(zlayers, (-Inf, z_top, _eps(geometry.stack.wafer.n)))
    push!(zlayers, (cbg_top, Inf, _eps(geometry.trench.n)))

    for j in 1:grid.Nz
        z_lo = grid.zp[j]; z_hi = grid.zp[j+1]; Δz = z_hi - z_lo
        base = zero(ComplexF64)
        for (z_bot, z_topL, ϵ) in zlayers
            lo = max(z_lo, z_bot); hi = min(z_hi, z_topL)
            hi > lo && (base += ((hi - lo) / Δz) * ϵ)
        end
        overlap = min(z_hi, cbg_top) - max(z_lo, cbg_bot)
        wband   = overlap > 0 ? overlap / Δz : 0.0
        if wband == 0
            @views ε[:, j] .= base
        else
            @views ε[:, j] .= base .+ wband .* ε_band
        end
    end

    return ε
end

function material_tensors(ϵ_grid, grid)
    ϵ_ρρ = zeros(ComplexF64, grid.Nρ, grid.Nz)
    ϵ_ρρ[:, 1]     .= ϵ_grid[:, 1]
    ϵ_ρρ[:, 2:end] .= 0.5 .* (ϵ_grid[:, 2:end] .+ ϵ_grid[:, 1:end-1])

    ϵ_zz = zeros(ComplexF64, grid.Nρ, grid.Nz)
    ϵ_zz[1, :]     .= ϵ_grid[1, :]
    ϵ_zz[2:end, :] .= 0.5 .* (ϵ_grid[2:end, :] .+ ϵ_grid[1:end-1, :])

    ϵ_ϕϕ = zeros(ComplexF64, grid.Nρ, grid.Nz)
    ϵ_ϕϕ[:, 1]     .= ϵ_zz[:, 1]
    ϵ_ϕϕ[:, 2:end] .= 0.5 .* (ϵ_zz[:, 2:end] .+ ϵ_zz[:, 1:end-1])

    Inv_Eρ = spdiagm(0 => 1.0 ./ vec(ϵ_ρρ))
    Inv_Eϕ = spdiagm(0 => 1.0 ./ vec(ϵ_ϕϕ))
    Inv_Ez = spdiagm(0 => 1.0 ./ vec(ϵ_zz))

    inv_ϵr = blockdiag(Inv_Eρ, Inv_Eϕ, Inv_Ez)
    inv_μr = spdiagm(0 => ones(ComplexF64, 3grid.N))
    return inv_ϵr, inv_μr
end