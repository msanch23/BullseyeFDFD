import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using BullseyeFDFD
using Printf, CSV, DataFrames

nClad = 1.00
nSiN  = 2.01066
nSiO2 = 1.45375
nSi   = 3.69476 - 1im*0.00482

rEye0  = 0.355000
tCBG0  = 0.17886
tSiO20 = 0.75057

λ0 = 0.780
NA = 0.4
Qfloor = 50

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

function eval_point(; rEye=rEye0, tCBG=tCBG0, tSiO2=tSiO20, λ_target=λ0)
    d = copy(design0); d[1] = rEye
    substrates = [(n = nSiO2, height = tSiO2),
                  (n = nSi,   height = Inf)]
    geometry = build_geometry(d, substrates, nClad, nSiN, tCBG)

    grid = conformal_grid(geometry, λ_target)
    ϵ = build_epsilon(geometry, grid)

    score(md) = md.CF - (md.Q ≤ Qfloor) * 10.0

    sim = solve_sim(grid, ϵ, λ_target, NA;
                    Nmodes=8, source=dipole(grid), select=score, report=false, plots=false)
    m = argmax(score, sim.modes)

    return (λ_nm = m.λ*1000, Q = m.Q, V = m.V_eff, Fp = m.Fp, CF = m.CF,
            Ceff_NA = m.Ceff_NA, Purcell = sim.driven.Purcell, η = m.η,
            gaussicity = m.gaussicity, Nρ = grid.Nρ, Nz = grid.Nz)
end

function run_sweep(param::Symbol, values)
    df = DataFrame()
    println("── sweeping $param ──")
    for v in sort(float.(collect(values)))
        r = eval_point(; Dict{Symbol,Any}(param => v)...)   # λ_target defaults to λ0
        @printf("  %s=%-7s → λ=%.2f Q=%.0f V=%.3f Fp=%.1f CF=%.3f Ceff=%.3f  %d×%d\n",
                param, string(v), r.λ_nm, r.Q, r.V, r.Fp, r.CF, r.Ceff_NA, r.Nρ, r.Nz)
        push!(df, (param=String(param), value=v,
                   λ_nm=r.λ_nm, Q=r.Q, V=r.V, Fp=r.Fp, CF=r.CF, Ceff_NA=r.Ceff_NA,
                   Purcell=r.Purcell, η=r.η, gaussicity=r.gaussicity, Nρ=r.Nρ, Nz=r.Nz))
    end
    return df
end

save(df, fname) = (CSV.write(joinpath(resultsdir, fname), df); println("Saved $fname\n"); nothing)

save(run_sweep(:rEye,  round.(0.310:0.010:0.400, digits=3)), "fdfd_rEye.csv")
save(run_sweep(:tCBG,  round.(0.160:0.005:0.210, digits=3)), "fdfd_tCBG.csv")
save(run_sweep(:tSiO2, round.(0.550:0.025:1.000, digits=3)), "fdfd_tSiO2.csv")
