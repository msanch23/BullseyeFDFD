struct CylindricalGrid
    R::Float64                  # Simulation Radial Size
    Z::Float64                  # Simulation Height Size 

    Nρ::Int                     # Number of Yee cells in the radial axis
    Nz::Int                     # Number of Yee cells in the vertical axis
    N::Int                      # Total number of Yee cells

    ρp::Vector{Float64}         # ρ primal locations    (length Nρ+1)
    zp::Vector{Float64}         # z primal locations    (length Nz+1)
    Δρp::Vector{Float64}        # Δρ primal widths      (length Nρ)
    Δzp::Vector{Float64}        # Δz primal widths      (length Nz)

    ρd::Vector{Float64}         # ρ dual locations      (length Nρ)
    zd::Vector{Float64}         # z dual locations      (length Nz)
    Δρd::Vector{Float64}        # Δρ dual widths        (length Nρ)
    Δzd::Vector{Float64}        # Δz dual widths        (length Nz)

    PMLρ::Float64               # μm thickness of PML cells in the radial direction
    PMLz_bot::Float64           # μm thickness of PML cells in the bottom vertical direction
    PMLz_top::Float64           # μm thickness of PML cells in the top vertical direction

    idx_PMLρ::Int               # Index of PML cells in the radial direction
    idx_PMLz_bot::Int           # Index of PML cells in the bottom vertical direction
    idx_PMLz_top::Int           # Index of PML cells in the top vertical direction

    cbg_bot::Float64            # μm location of the bottom of the Bullseye cavity
    cbg_cen::Float64            # μm location of the center of the Bullseye cavity
    cbg_top::Float64            # μm location of the top of the Bullseye cavity

    idx_cbg_bot::Int            # Index of the bottom of the Bullseye cavity
    idx_cbg_cen::Int            # Index of the center of the Bullseye cavity
    idx_cbg_top::Int            # Index of the top of the Bullseye cavity

    total_cbg_radius::Float64   # Total Bullseye radius (includes buffer)
    eye_radius::Float64         # Radius of the central eye feature
end;

function suggest_Δ(geometry, λ, λPts, geoPts)
    optΔ(n) = λ / (real(n) * λPts)

    min_Δρ = Inf
    haskey(geometry, :slot) && (min_Δρ = min(min_Δρ, optΔ(geometry.slot.n), geometry.slot.width / geoPts))
    for feat in (geometry.eye, geometry.buffer)
        min_Δρ = min(min_Δρ, optΔ(feat.n), feat.width / geoPts)
    end

    for grp in (geometry.trench, geometry.ring)
        for w in grp.width
            min_Δρ = min(min_Δρ, optΔ(grp.n), w / geoPts)
        end
    end

    min_Δz = min(optΔ(geometry.ring.n), geometry.tCBG / geoPts)
    for sub in Base.front(geometry.stack)
        min_Δz = min(min_Δz, optΔ(sub.n), sub.height / geoPts)
    end

    return min_Δρ, min_Δz
end;

function stretch_grid(L, Δ, stretch_ratio, max_stretch_Δ)
    @assert Δ > 0 && max_stretch_Δ > 0 "stretch_grid needs positive Δ and max_stretch_Δ"

    stretched_Δ = Float64[]
    current_L = 0.0
    current_Δ = Δ

    while current_L < L
        if current_L + current_Δ > L
            remaining_L = L - current_L

            if remaining_L < 0.5 * current_Δ && !isempty(stretched_Δ)
                stretched_Δ[end] += remaining_L
            elseif remaining_L > 0
                push!(stretched_Δ, remaining_L)
            end
            break
        end

        push!(stretched_Δ, current_Δ)
        current_L += current_Δ
        current_Δ = min(current_Δ * stretch_ratio, max_stretch_Δ)
    end

    return stretched_Δ
end;

function conformal_grid(geometry, λ;
                        λPts=20, geoPts=16,
                        PMLρ=0.5, PMLz_bot=0.5, PMLz_top=0.5,
                        padρ = 0.5, padz_bot = 0.5, padz_top = 0.5,
                        stretch_ratio_Δρ = 1.15, stretch_ratio_Δz_bot = 1.15, stretch_ratio_Δz_top = 1.15,
                        max_stretch_Δρ = 0.0, max_stretch_Δz_bot = 0.0, max_stretch_Δz_top = 0.0)

    target_Δρ, target_Δz = suggest_Δ(geometry, λ, λPts, geoPts)
    max_Δρ_limit = iszero(max_stretch_Δρ) ? 3*target_Δρ : max_stretch_Δρ
    max_Δz_bot_limit = iszero(max_stretch_Δz_bot) ? 3*target_Δz : max_stretch_Δz_bot
    max_Δz_top_limit = iszero(max_stretch_Δz_top) ? 3*target_Δz : max_stretch_Δz_top

    vector_Δρ = Float64[]
    pixelate!(v, w) = append!(v, fill(w / max(1, ceil(Int, w/target_Δρ)),
                                          max(1, ceil(Int, w/target_Δρ))))
    haskey(geometry, :slot) && pixelate!(vector_Δρ, geometry.slot.width)
    pixelate!(vector_Δρ, geometry.eye.width)
    for i in eachindex(geometry.ring.width)
        pixelate!(vector_Δρ, geometry.trench.width[i])
        pixelate!(vector_Δρ, geometry.ring.width[i])
    end
    pixelate!(vector_Δρ, geometry.buffer.width)
    slot_w = haskey(geometry, :slot) ? geometry.slot.width : 0.0
    total_cbg_radius = slot_w + geometry.eye.width + sum(geometry.trench.width) +
                       sum(geometry.ring.width) + geometry.buffer.width
    eye_radius = geometry.eye.width

    remaining_ρ = (total_cbg_radius + padρ + PMLρ) - sum(vector_Δρ)
    if remaining_ρ > 0
        start_Δ = isempty(vector_Δρ) ? target_Δρ : vector_Δρ[end]
        append!(vector_Δρ, stretch_grid(remaining_ρ, start_Δ, stretch_ratio_Δρ, max_Δρ_limit))
    end
    Nρ = length(vector_Δρ)
    ρp = vcat(0.0, cumsum(vector_Δρ))
    R = ρp[end]
    ρd = ρp[1:end-1] .+ 0.5 .* vector_Δρ
    Δρd = vcat(vector_Δρ[1], ρd[2:end] - ρd[1:end-1])

    dev_height = geometry.tCBG
    N_dev = max(1, ceil(Int, dev_height / target_Δz))
    vector_bullseye_Δz = fill(dev_height / N_dev, N_dev)

    vector_bot_Δz = Float64[]
    finite_substrates = Base.front(geometry.stack) 

    for sub in finite_substrates
        sub_Δz_target = λ / (real(sub.n) * λPts)
        start_Δ = isempty(vector_bot_Δz) ? vector_bullseye_Δz[end] : vector_bot_Δz[1]
        cells = stretch_grid(sub.height, start_Δ, stretch_ratio_Δz_bot, sub_Δz_target)
        prepend!(vector_bot_Δz, reverse(cells))
    end

    L_bot_remaining = padz_bot + PMLz_bot
    if L_bot_remaining > 0
        start_Δ = isempty(vector_bot_Δz) ? target_Δz : vector_bot_Δz[1]
        vec = stretch_grid(L_bot_remaining, start_Δ, stretch_ratio_Δz_bot, max_Δz_bot_limit)
        prepend!(vector_bot_Δz, reverse(vec))
    end

    vector_top_Δz = Float64[]
    L_top_total = padz_top + PMLz_top
    if L_top_total > 0
        start_Δ = isempty(vector_bullseye_Δz) ? target_Δz : vector_bullseye_Δz[end]
        vector_top_Δz = stretch_grid(L_top_total, start_Δ, stretch_ratio_Δz_top, max_Δz_top_limit)
    end

    vector_Δz = vcat(vector_bot_Δz, vector_bullseye_Δz, vector_top_Δz)
    Nz = length(vector_Δz)
    zp = vcat(0.0, cumsum(vector_Δz))
    Z = zp[end]
    zd = zp[1:end-1] .+ 0.5 .* vector_Δz
    Δzd = vcat(zd[1], zd[2:end] - zd[1:end-1])
    N = Nρ * Nz

    cbg_bot = sum(vector_bot_Δz)
    cbg_top = cbg_bot + sum(vector_bullseye_Δz)
    cbg_cen = (cbg_bot + cbg_top) / 2.0
    idx_cbg_bot = length(vector_bot_Δz) + 1
    idx_cbg_top = idx_cbg_bot + length(vector_bullseye_Δz) - 1
    idx_cbg_cen = idx_cbg_bot + div(length(vector_bullseye_Δz), 2)

    idx_PMLρ = argmin(abs.(ρd .- (R - PMLρ)))
    idx_PMLz_bot = argmin(abs.(zd .- PMLz_bot))
    idx_PMLz_top = argmin(abs.(zd .- (Z - PMLz_top)))

    return CylindricalGrid(R, Z, Nρ, Nz, N,
                           ρp, zp, vector_Δρ, vector_Δz,
                           ρd, zd, Δρd, Δzd,
                           PMLρ, PMLz_bot, PMLz_top,
                           idx_PMLρ, idx_PMLz_bot, idx_PMLz_top,
                           cbg_bot, cbg_cen, cbg_top,
                           idx_cbg_bot, idx_cbg_cen, idx_cbg_top,
                           total_cbg_radius, eye_radius)
end;