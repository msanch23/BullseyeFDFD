module BullseyeFDFD

using LinearAlgebra
using SparseArrays
using Plots
using Arpack
using Printf
using SpecialFunctions
using Optim
using Dates
using LineSearches
using LBFGSB: L_BFGS_B, setulb

include("geometry.jl")     # build_geometry: design vector → role blocks (eye/trench/ring/buffer) + stack
include("grid.jl")         # CylindricalGrid + conformal_grid ((ρ, z) mesh: feature edges on nodes, stretched padding)
include("operators.jl")    # s_vectors (PML stretching) + CeChOperators (curls, ρ=0 axis, m)
include("materials.jl")    # build_epsilon (Yee-voxel ε, ramped walls) + material_tensors
include("solver.jl")       # solve_sim: shift-invert eigensolve + driven dipole, mode selection
include("postprocess.jl")  # mode/driven analysis, V/Purcell/CF, flux, far-field, plotting
include("adjoint.jl")      # eigenvalue/eigenvector adjoints, FOM terms, adj_opt driver

export build_geometry
export conformal_grid
export build_epsilon
export show_sim
export CeChOperators
export solve_sim, dipole
export analyze_mode, report_modes, confinement_factor, analyze_driven, report_driven
export plot_efield, plot_farfield, plot_kspace, plot_farfield_kspace
export calc_farfield, calc_flux, calc_H_fields
export adj_opt, λ_target, Q_target, V_target, Fp_target, ηG_target, η_target, FOMTerm
export BackTracking, HagerZhang

end
