module BullseyeFDFD

using LinearAlgebra
using SparseArrays
using Plots
using Arpack
using Printf
using SpecialFunctions
using Optim
using LineSearches

include("geometry.jl")     # Feature/Layer types + build_geometry, design ↔ params bridge
include("grid.jl")         # cylindrical (r, z) mesh + coordinate stretching / PML profiles
include("operators.jl")    # sparse derivative / curl operators on the grid
include("materials.jl")    # ε, μ, σ assignment onto the grid
include("solver.jl")       # system assembly + linear solve
include("postprocess.jl")  # field extraction, modes, fluxes, plotting helpers
include("adjoint.jl")      # eigenvalue-sensitivity adjoint + FOM terms + optimizer driver

export build_geometry
export conformal_grid
export build_epsilon
export show_sim          
export CeChOperators
export solve_sim, dipole
export analyze_mode, report_modes, confinement_factor, analyze_driven, report_driven
export plot_efield, plot_farfield, plot_kspace, plot_farfield_kspace
export calc_farfield, calc_flux, calc_H_fields
export adj_opt, λ_target, Q_target, FOMTerm
export BackTracking, HagerZhang

end
