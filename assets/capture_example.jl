# Re-run the README walkthrough (examples/sanchezBullseye.jl) headless, unchanged, and save what the
# README shows: the report text → assets/readme/*.txt, the figures → assets/*.png.
# Then `julia assets/update_readme.jl` copies the text into README.md.
#
#   julia --project=. assets/capture_example.jl
ENV["GKSwstype"] = "100"
using BullseyeFDFD, Dates
using Plots: savefig

const ROOT   = dirname(@__DIR__)
const OUT    = joinpath(ROOT, "assets", "readme")
const ASSETS = joinpath(ROOT, "assets")
mkpath(OUT)

# run the example with its stdout captured; its globals (grid, ϵ, λ0, NA, sim) stay in Main
logfile = tempname()
open(logfile, "w") do io
    redirect_stdout(() -> include(joinpath(ROOT, "examples", "sanchezBullseye.jl")), io)
end
lines = split(read(logfile, String), '\n')

# report_modes prints a title and a table between two "====" rules; report_driven prints two lines
tables, drivens = String[], String[]
let i = 1
    while i <= length(lines)
        if strip(lines[i]) == "Eigenmode Analysis Report"
            j = findnext(l -> startswith(l, "===="), lines, findnext(l -> startswith(l, "===="), lines, i) + 1)
            push!(tables, join(lines[i:j], '\n')); i = j + 1
        elseif startswith(lines[i], "Purcell factor = ")
            push!(drivens, join(lines[i:i+1], '\n')); i += 2
        else
            i += 1
        end
    end
end
length(tables) == 3 && length(drivens) == 2 ||
    error("expected 3 eigenmode reports and 2 driven reports from the example, got $(length(tables)) and $(length(drivens))")

repo = pkgdir(BullseyeFDFD)
prov = string("# provenance: git=", readchomp(`git -C $repo rev-parse HEAD`),
              " src_clean=", success(`git -C $repo diff --quiet HEAD -- src Project.toml`),
              " version=", pkgversion(BullseyeFDFD), " grid=", grid.Nρ, "×", grid.Nz,
              " (conformal_grid defaults) date=", today())
save(name, text) = write(joinpath(OUT, name), prov, '\n', text, '\n')

# README order: Nmodes=3; + NA; driven at λ_drive; NA + Nmodes=3 + dipole at the mode
save("sanchez_eig.txt",        tables[1])
save("sanchez_eig_na.txt",     tables[2])
save("sanchez_driven.txt",     drivens[1])
save("sanchez_eig_driven.txt", tables[3] * "\n\n" * drivens[2])

# figures: the cross-section, and the plots solve_sim draws for its last call (mode chosen as solve_sim does)
k  = BullseyeFDFD._select_index(sim.modes, λ0, nothing, grid, sim.inv_ϵr)
md = sim.modes[k]
savefig(show_sim(ϵ, grid),                                 joinpath(ASSETS, "crosssection.png"))
savefig(plot_efield(grid, md; log_scale=false, label="Mode $k"), joinpath(ASSETS, "eigenmode.png"))
savefig(plot_farfield_kspace(md.θs, md.farfield, NA),      joinpath(ASSETS, "farfield_kspace.png"))
savefig(plot_efield(grid, sim.driven; log_scale=true,
                    label="Driven @ $(round(sim.driven.λ*1000; digits=1)) nm"), joinpath(ASSETS, "driven_efield.png"))
println("saved 4 report blocks → $OUT, 4 figures → $ASSETS")
