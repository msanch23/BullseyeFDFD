_eps(n) = ComplexF64(n)^2

# Share of the voxel [a, b] lying below x (ρ: by area, z: by length) and its antiderivative in x.
# _below ramps each step over a width w, which makes every share C¹ in the interface position.
_Gρ(x, a, b) = x <= a ? 0.0 : x >= b ? 1.0 : (x^2 - a^2) / (b^2 - a^2)
_Gz(x, a, b) = x <= a ? 0.0 : x >= b ? 1.0 : (x - a) / (b - a)
_Pρ(x, a, b) = x <= a ? 0.0 : x >= b ? (b^3/3 - a^2*b + 2a^3/3) / (b^2 - a^2) + (x - b) :
                                       (x^3/3 - a^2*x + 2a^3/3) / (b^2 - a^2)
_Pz(x, a, b) = x <= a ? 0.0 : x >= b ? (b - a)/2 + (x - b) : (x - a)^2 / (2(b - a))

_below(G, P, x, a, b, w) = w == 0 ? G(x, a, b) : (P(x + w/2, a, b) - P(x - w/2, a, b)) / w

function _fractions(G, P, (a, b), interfaces, w)
    f = Vector{Float64}(undef, length(interfaces) + 1); prev = 0.0
    for k in eachindex(interfaces)
        g = _below(G, P, interfaces[k], a, b, w); f[k] = g - prev; prev = g
    end
    f[end] = 1 - prev
    return f
end

# radial order slot→eye→(trench[i], ring[i])→buffer, then bulk (eye.n) out to R
function _radial_profile(geometry)
    has_slot = haskey(geometry, :slot)
    widths   = has_slot ? Float64[geometry.slot.width, geometry.eye.width]        : Float64[geometry.eye.width]
    eps_seg  = has_slot ? ComplexF64[_eps(geometry.slot.n), _eps(geometry.eye.n)] : ComplexF64[_eps(geometry.eye.n)]
    for i in eachindex(geometry.ring.width)
        push!(widths, geometry.trench.width[i]); push!(eps_seg, _eps(geometry.trench.n))
        push!(widths, geometry.ring.width[i]);   push!(eps_seg, _eps(geometry.ring.n))
    end
    push!(widths, geometry.buffer.width); push!(eps_seg, _eps(geometry.buffer.n))
    push!(eps_seg, _eps(geometry.eye.n))
    return cumsum(widths), eps_seg
end

# layers bottom→top: wafer, substrates (stacked down from grid.cbg_bot), band, cladding.
# The band's ε varies with ρ, so its entry in eps_lay is only a placeholder.
function _vertical_profile(geometry, grid)
    z_top = grid.cbg_bot; tops = Float64[]; eps_sub = ComplexF64[]
    for sub in Base.front(geometry.stack)
        push!(tops, z_top); push!(eps_sub, _eps(sub.n)); z_top -= sub.height
    end
    z_interfaces = vcat(z_top, reverse(tops), grid.cbg_bot + geometry.tCBG)
    eps_lay      = vcat(_eps(geometry.stack.wafer.n), reverse(eps_sub), NaN, _eps(geometry.trench.n))
    return z_interfaces, eps_lay, length(eps_sub) + 2
end

# band ε over each ρ-voxel: arithmetic, or harmonic for the component normal to the ρ-walls
_band(voxels, ρ_interfaces, eps_seg, w; harmonic=false) = map(voxels) do v
    f = _fractions(_Gρ, _Pρ, v, ρ_interfaces, w)
    harmonic ? 1 / sum(f ./ eps_seg) : sum(f .* eps_seg)
end

# per z-voxel: (Σ over the non-band layers of f·ε, or of f/ε for the component normal to the
# z-walls; the band's fraction)
_layers(voxels, z_interfaces, eps_lay, iband, w; inverse=false) = map(voxels) do v
    f = _fractions(_Gz, _Pz, v, z_interfaces, w); acc = zero(ComplexF64)
    for k in eachindex(f)
        (k == iband || f[k] == 0) && continue
        acc += f[k] * (inverse ? 1 / eps_lay[k] : eps_lay[k])
    end
    (acc, f[iband])
end

_compose(lay, band) = [lay[j][1] + lay[j][2] * band[i] for i in eachindex(band), j in eachindex(lay)]

# Ramp widths: the smallest cell an interface can occupy, so each voxel's kernel is a hat (w = Δ)
# or a trapezoid (w < Δ) and never reaches past the neighbouring voxel. Depends on the grid only,
# so it stays fixed while the design moves on a frozen grid.
_default_ramp(grid) = (minimum(grid.Δρp[grid.ρd .< grid.total_cbg_radius]),
                       minimum(grid.Δzp[grid.PMLz_bot .< grid.zd .< grid.cbg_top]))

# cell: sharp cell-averaged map (geometry, plotting). ρρ, ϕϕ, zz: ε at the Eρ, Eϕ, Ez Yee nodes,
# each averaged over that component's own voxel with ramped walls, harmonic across a wall for the
# normal component and arithmetic for the tangential ones.
function build_epsilon(geometry, grid; ramp = _default_ramp(grid))
    wρ, wz = ramp
    ρ_interfaces, eps_seg        = _radial_profile(geometry)
    z_interfaces, eps_lay, iband = _vertical_profile(geometry, grid)

    ρ_cell = [(grid.ρp[i], grid.ρp[i+1]) for i in 1:grid.Nρ]
    ρ_node = [(i == 1 ? 0.0 : grid.ρd[i-1], grid.ρd[i]) for i in 1:grid.Nρ]
    z_cell = [(grid.zp[j], grid.zp[j+1]) for j in 1:grid.Nz]
    z_node = [(j == 1 ? grid.zp[1] : grid.zd[j-1], grid.zd[j]) for j in 1:grid.Nz]

    band(v, w; kw...)   = _band(v, ρ_interfaces, eps_seg, w; kw...)
    layers(v, w; kw...) = _layers(v, z_interfaces, eps_lay, iband, w; kw...)
    lay_node = layers(z_node, wz)

    cell = _compose(layers(z_cell, 0.0), band(ρ_cell, 0.0))
    ρρ   = _compose(lay_node, band(ρ_cell, wρ; harmonic=true))
    ϕϕ   = _compose(lay_node, band(ρ_node, wρ))
    zz   = 1 ./ _compose(layers(z_cell, wz; inverse=true), 1 ./ band(ρ_node, wρ))
    return (; cell, ρρ, ϕϕ, zz)
end

_cell(ϵ::AbstractMatrix) = ϵ
_cell(ϵ::NamedTuple)     = ϵ.cell

_yee_eps(ϵ::NamedTuple, grid) = vcat(vec(ϵ.ρρ), vec(ϵ.ϕϕ), vec(ϵ.zz))

function _yee_eps(ϵ_grid::AbstractMatrix, grid)
    ϵ_ρρ = zeros(ComplexF64, grid.Nρ, grid.Nz)
    ϵ_ρρ[:, 1]     .= ϵ_grid[:, 1]
    ϵ_ρρ[:, 2:end] .= 0.5 .* (ϵ_grid[:, 2:end] .+ ϵ_grid[:, 1:end-1])

    ϵ_zz = zeros(ComplexF64, grid.Nρ, grid.Nz)
    ϵ_zz[1, :]     .= ϵ_grid[1, :]
    ϵ_zz[2:end, :] .= 0.5 .* (ϵ_grid[2:end, :] .+ ϵ_grid[1:end-1, :])

    ϵ_ϕϕ = zeros(ComplexF64, grid.Nρ, grid.Nz)
    ϵ_ϕϕ[:, 1]     .= ϵ_zz[:, 1]
    ϵ_ϕϕ[:, 2:end] .= 0.5 .* (ϵ_zz[:, 2:end] .+ ϵ_zz[:, 1:end-1])

    return vcat(vec(ϵ_ρρ), vec(ϵ_ϕϕ), vec(ϵ_zz))
end

function material_tensors(ϵ, grid)
    inv_ϵr = spdiagm(0 => 1.0 ./ _yee_eps(ϵ, grid))
    inv_μr = spdiagm(0 => ones(ComplexF64, 3grid.N))
    return inv_ϵr, inv_μr
end