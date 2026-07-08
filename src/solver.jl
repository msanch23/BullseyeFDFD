function _check_source(source, m)
    need = source.pol === :z ? 0 : 1
    m == need || error("pol=$(source.pol) requires m=$need, got m=$m")
end

function _unpack_field(e_field, grid)
    N = grid.N
    Eρ = reshape(e_field[1:N],        grid.Nρ, grid.Nz)
    Eϕ = reshape(e_field[N+1:2N],     grid.Nρ, grid.Nz)
    Ez = reshape(e_field[2N+1:3N],    grid.Nρ, grid.Nz)
    return Eρ, Eϕ, Ez
end

function _assemble(grid, ϵ, λ; m=1, pml_kwargs...)
    Ce, Ch        = CeChOperators(grid, λ; m=m, pml_kwargs...)
    inv_ϵr, inv_μr = material_tensors(ϵ, grid)
    A = inv_ϵr * Ch * inv_μr * Ce
    return A, inv_ϵr, Ce 
end

function _solve_eig(A, grid, λ_target, Nmodes, tol, maxiter)
    k0_sq = (2π / λ_target)^2 + 0im
    vals, vecs = eigs(A; nev=Nmodes, sigma=k0_sq, which=:LM,
                      tol=tol, maxiter=maxiter)

    modes = Vector{NamedTuple}(undef, Nmodes)
    for idx in 1:Nmodes
        k_sq = vals[idx]
        k    = sqrt(k_sq)
        λ_res = 2π / real(k)
        Q     = real(k) / (2 * abs(imag(k)))

        e_field    = vecs[:, idx]
        Eρ, Eϕ, Ez = _unpack_field(e_field, grid)
        modes[idx] = (; λ=λ_res, Q, k_sq, raw_E=e_field, Eρ, Eϕ, Ez)
    end
    return modes
end

function _solve_driven(A, inv_ϵr, grid, λ, source)
    k0    = 2π / λ
    k0_sq = k0^2 + 0im

    L = A - k0_sq*I
    F = lu(L)

    block  = source.pol === :ρ ? 0 : source.pol === :ϕ ? grid.N : 2grid.N
    flat   = block + (source.z_idx - 1) * grid.Nρ + source.ρ_idx
    J = zeros(ComplexF64, 3grid.N)
    J[flat] = 1.0
    b = 1im * k0 * (inv_ϵr * J)

    e_field    = F \ b
    Eρ, Eϕ, Ez = _unpack_field(e_field, grid)

    E_squared = abs2.(Eρ) .+ abs2.(Eϕ) .+ abs2.(Ez)
    normE     = sqrt.(E_squared)

    return (; λ, k_sq=k0_sq, raw_E=e_field, Eρ, Eϕ, Ez,
              E_squared, normE, J, source, F)
end

function solve_sim(grid, ϵ, λ, NA=nothing;
                   Nmodes=nothing, source=nothing, m=1,
                   select=nothing, report=true, purcell=true, plots = true,
                   tol=1e-6, maxiter=3000, pml_kwargs...)
    @assert Nmodes !== nothing || source !== nothing "provide Nmodes and/or source"
    source === nothing || _check_source(source, m)

    A, inv_ϵr, Ce = _assemble(grid, ϵ, λ; m=m, pml_kwargs...)

    modes = nothing
    if Nmodes !== nothing
        raw   = _solve_eig(A, grid, λ, Nmodes, tol, maxiter)
        modes = [analyze_mode(md, grid, ϵ, inv_ϵr) for md in raw]
        if NA !== nothing
            mon   = default_monitors(grid)
            modes = [(; md..., analyze_collection(md, grid, Ce; NA, monitors=mon)...)
                     for md in modes]
        end
        report && report_modes(modes)
    end


    score = select === nothing ? (md -> md.CF) : select

    driven = nothing
    if source !== nothing
        drive_mode = modes === nothing ? nothing : argmax(score, modes)
        λ_d = source.λ   !== nothing ? source.λ    :
              drive_mode !== nothing ? drive_mode.λ : λ

        A_d, iϵ_d = λ_d == λ ? (A, inv_ϵr) : _assemble(grid, ϵ, λ_d; m=m, pml_kwargs...)
        driven = analyze_driven(_solve_driven(A_d, iϵ_d, grid, λ_d, source), grid)

        if purcell
            ϵ_bulk    = fill(ϵ[source.ρ_idx, source.z_idx], grid.Nρ, grid.Nz)
            A_b, iϵ_b = _assemble(grid, ϵ_bulk, λ_d; m=m, pml_kwargs...)
            d_bulk    = _solve_driven(A_b, iϵ_b, grid, λ_d, source)
            driven    = (; driven..., Purcell = source_power(grid, driven) /
                                                source_power(grid, d_bulk))
        end

        report && report_driven(driven)
    end

    if plots
        if modes !== nothing
            k = argmax(i -> score(modes[i]), eachindex(modes))
            plot_efield(grid, modes[k]; log_scale=false, label="Mode $k")
            NA !== nothing && plot_farfield_kspace(modes[k].θs, modes[k].farfield, NA)
        end
        driven !== nothing && plot_efield(grid, driven; log_scale=true,
                                          label="Driven @ $(round(driven.λ*1000; digits=1)) nm")
    end

    return (; modes, driven, A, inv_ϵr, Ce)
end

dipole(grid; λ=nothing, ρ_idx=2, z_idx=grid.idx_cbg_cen, pol=:ρ) =
    (; λ, ρ_idx, z_idx, pol)
dipole(grid, λ; kwargs...) = dipole(grid; λ=float(λ), kwargs...)
