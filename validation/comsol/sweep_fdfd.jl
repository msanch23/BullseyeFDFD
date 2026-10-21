import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using BullseyeFDFD
using Printf, CSV, DataFrames, Dates

nClad = 1.00
nSiN  = 2.01066
nSiO2 = 1.45375
nSi   = 3.69476 - 1im*0.00482

rEye0  = 0.355000
tCBG0  = 0.17886
tSiO20 = 0.75057

λ0 = 0.780
NA = 0.4
gridkw = (; λPts=36, geoPts=24)      # the accurate setting (validation/convergence, sections 8–9); recorded in every row

design0 = [rEye0,               # eye  (index 1 — swept via rEye)
           0.100000, 0.370000,  # trench 1, ring 1
           0.100000, 0.370000,  # trench 2, ring 2
           0.100000, 0.370000,  # trench 3, ring 3
           0.100000, 0.370000,  # trench 4, ring 4
           0.100000, 0.370000,  # trench 5, ring 5
           0.100000, 0.370000,  # trench 6, ring 6
           0.100000, 0.370000,  # trench 7, ring 7
           0.100000, 0.370000,  # trench 8, ring 8
           0.100000, 0.370000,  # trench 9, ring 9
           0.100000, 0.370000,  # trench 10, ring 10
           0.100000]            # buffer

resultsdir = joinpath(@__DIR__, "results")
mkpath(resultsdir)

# provenance stamped on every row: the commit, whether src/ + Project.toml match it, version, date
repo = pkgdir(BullseyeFDFD)
prov = (git       = readchomp(`git -C $repo rev-parse HEAD`),
        src_clean = success(`git -C $repo diff --quiet HEAD -- src Project.toml`),
        version   = string(pkgversion(BullseyeFDFD)),
        date      = string(today()))

function eval_point(; rEye=rEye0, tCBG=tCBG0, tSiO2=tSiO20, λ_target=λ0)
    d = copy(design0); d[1] = rEye
    substrates = [(n = nSiO2, height = tSiO2),
                  (n = nSi,   height = Inf)]
    geometry = build_geometry(d, substrates, nClad, nSiN, tCBG)

    grid = conformal_grid(geometry, λ_target; gridkw...)
    ϵ = build_epsilon(geometry, grid)

    # the package's default mode choice (Q ≥ 30, PML fraction ≤ 0.5, Q weighted by distance from λ_target);
    # solve_sim drives the dipole at that same mode
    sim = solve_sim(grid, ϵ, λ_target, NA;
                    Nmodes=8, source=dipole(grid), report=false, plots=false)
    m = sim.modes[BullseyeFDFD._select_index(sim.modes, λ_target, nothing, grid, sim.inv_ϵr)]

    return (λ_nm = m.λ*1000, Q = m.Q, V = m.V_eff, Fp = m.Fp, CF = m.CF,
            Ceff_NA = m.Ceff_NA, Purcell = sim.driven.Purcell, η = m.η,
            gaussicity = m.gaussicity, pml_frac = m.pml_frac, Nρ = grid.Nρ, Nz = grid.Nz)
end

function run_sweep(param::Symbol, values)
    df = DataFrame()
    println("── sweeping $param ──")
    for v in sort(float.(collect(values)))
        r = eval_point(; Dict{Symbol,Any}(param => v)...)   # λ_target defaults to λ0
        @printf("  %s=%-7s → λ=%.2f Q=%.0f V=%.3f Fp=%.1f CF=%.3f Ceff=%.3f PML=%.3f  %d×%d\n",
                param, string(v), r.λ_nm, r.Q, r.V, r.Fp, r.CF, r.Ceff_NA, r.pml_frac, r.Nρ, r.Nz)
        push!(df, (param=String(param), value=v,
                   λ_nm=r.λ_nm, Q=r.Q, V=r.V, Fp=r.Fp, CF=r.CF, Ceff_NA=r.Ceff_NA,
                   Purcell=r.Purcell, η=r.η, gaussicity=r.gaussicity, Nρ=r.Nρ, Nz=r.Nz,
                   pml_frac=r.pml_frac, gridkw..., prov...))
    end
    return df
end

save(df, fname) = (CSV.write(joinpath(resultsdir, fname), df); println("Saved $fname\n"); nothing)

# `julia sweep_fdfd.jl tSiO2` re-runs one sweep; no arguments runs all three
sweeps = Dict(:rEye  => round.(0.310:0.010:0.400, digits=3),
              :tCBG  => round.(0.160:0.005:0.210, digits=3),
              :tSiO2 => round.(0.550:0.025:1.000, digits=3))
for p in (isempty(ARGS) ? [:rEye, :tCBG, :tSiO2] : Symbol.(ARGS))
    save(run_sweep(p, sweeps[p]), "fdfd_$p.csv")
end
