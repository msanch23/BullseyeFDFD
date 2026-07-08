function s_vectors(grid::CylindricalGrid, λ; p=3.5, R_target=1e-12)
    calc_σ_max = (L) -> -(p + 1) * log(R_target) / (2 * L)

    function calc_S(u, L, σ_max)
        if L <= 0 || u <= 0; return 1.0 + 0im; end
        d = u / L
        σ = σ_max * d^p
        return 1.0 - (σ * λ / (2π)) * im
    end

    Sρp = ones(ComplexF64, grid.Nρ)
    Sρd = ones(ComplexF64, grid.Nρ)

    if grid.PMLρ > 0
        σ_max_ρ = calc_σ_max(grid.PMLρ)
        ρ_start = grid.R - grid.PMLρ
        for i in 1:grid.Nρ
            if grid.ρp[i] > ρ_start
                Sρp[i] = calc_S(grid.ρp[i] - ρ_start, grid.PMLρ, σ_max_ρ)
            end
            if grid.ρd[i] > ρ_start
                Sρd[i] = calc_S(grid.ρd[i] - ρ_start, grid.PMLρ, σ_max_ρ)
            end
        end
    end

    Szp = ones(ComplexF64, grid.Nz)
    Szd = ones(ComplexF64, grid.Nz)

    if grid.PMLz_bot > 0
        σ_max_bot = calc_σ_max(grid.PMLz_bot)
        z_end = grid.PMLz_bot
        for i in 1:grid.Nz
            if grid.zp[i] < z_end
                Szp[i] = calc_S(z_end - grid.zp[i], grid.PMLz_bot, σ_max_bot)
            end
            if grid.zd[i] < z_end
                Szd[i] = calc_S(z_end - grid.zd[i], grid.PMLz_bot, σ_max_bot)
            end
        end
    end

    if grid.PMLz_top > 0
        σ_max_top = calc_σ_max(grid.PMLz_top)
        z_start = grid.Z - grid.PMLz_top
        for i in 1:grid.Nz
            if grid.zp[i] > z_start
                Szp[i] = calc_S(grid.zp[i] - z_start, grid.PMLz_top, σ_max_top)
            end
            if grid.zd[i] > z_start
                Szd[i] = calc_S(grid.zd[i] - z_start, grid.PMLz_top, σ_max_top)
            end
        end
    end

    return Sρp, Sρd, Szp, Szd
end;

function CeChOperators(grid::CylindricalGrid, λ; m=1, pml_kwargs...)
    Sρp, Sρd, Szp, Szd = s_vectors(grid, λ; pml_kwargs...)

    inv_Δzp = 1.0 ./ (grid.Δzp .* Szd)
    Dz_E = spdiagm(0 => -inv_Δzp, 1 => inv_Δzp[1:end-1])

    inv_Δzd = 1.0 ./ (grid.Δzd .* Szp)
    Dz_H = spdiagm(0 => inv_Δzd, -1 => -inv_Δzd[2:end])

    inv_Δρp = 1.0 ./ (grid.Δρp .* Sρd)
    Dρ_E = spdiagm(0 => -inv_Δρp, 1 => inv_Δρp[1:end-1])

    inv_Δρd = 1.0 ./ (grid.Δρd .* Sρp)
    d0 = copy(inv_Δρd)
    d0[1] *= 2.0
    Dρ_H = spdiagm(0 => d0, -1 => -inv_Δρd[2:end])

    Iρ = sparse(I, grid.Nρ, grid.Nρ)
    Iz = sparse(I, grid.Nz, grid.Nz)

    Dz_E_2D = kron(Dz_E, Iρ)
    Dz_H_2D = kron(Dz_H, Iρ)
    Dρ_E_2D = kron(Iz, Dρ_E)
    Dρ_H_2D = kron(Iz, Dρ_H)

    ρp_vec = grid.ρp[1:end-1]
    ρd_vec = grid.ρd
    inv_ρp_vec = [r == 0 ? 0.0 : 1.0/r for r in ρp_vec]
    inv_ρd_vec = 1.0 ./ ρd_vec

    ρp_2D     = kron(Iz, spdiagm(0 => ρp_vec))
    ρd_2D     = kron(Iz, spdiagm(0 => ρd_vec))
    inv_ρp_2D = kron(Iz, spdiagm(0 => inv_ρp_vec))
    inv_ρd_2D = kron(Iz, spdiagm(0 => inv_ρd_vec))

    Dρ_ρE = inv_ρd_2D * Dρ_E_2D * ρp_2D
    Dρ_ρH = inv_ρp_2D * Dρ_H_2D * ρd_2D

    Z = spzeros(ComplexF64, grid.N, grid.N)

    Ce = [ Z                -Dz_E_2D  -1im*m*inv_ρp_2D ;
           Dz_E_2D           Z        -Dρ_E_2D         ;
           1im*m*inv_ρd_2D   Dρ_ρE     Z               ]

    Ch = [ Z                -Dz_H_2D  -1im*m*inv_ρd_2D ;
           Dz_H_2D           Z        -Dρ_H_2D         ;
           1im*m*inv_ρp_2D   Dρ_ρH     Z               ]

    return Ce, Ch
end;
