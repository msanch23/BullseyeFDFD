# Theory notes

Derivations and background for the cylindrical FDFD solver, versioned alongside
the code they describe. One note per topic, each kept close to the `src/` file it
explains. Symbols are defined once in [`notation.md`](notation.md); every other
note links there instead of redefining them.

## Reading order

1. [`formulation.md`](formulation.md)
2. [`geometry.md`](geometry.md)
3. [`discretization.md`](discretization.md)
4. [`pml.md`](pml.md)
5. [`materials.md`](materials.md)
6. [`system-assembly.md`](system-assembly.md)
7. [`postprocessing.md`](postprocessing.md)
8. [`adjoint.md`](adjoint.md)
9. [`validation.md`](validation.md)

## Notes

| Note | Covers | Code | Status |
|---|---|---|---|
| [`formulation.md`](formulation.md) | Maxwell's equations in (r, φ, z), azimuthal mode decomposition (e^{imφ}), and the reduction to a 2D (r, z) problem | — | drafted |
| [`geometry.md`](geometry.md) | bullseye geometry parametrization and role blocks | `src/geometry.jl` | planned |
| [`discretization.md`](discretization.md) | Yee-cell layout in cylindrical coordinates, the discrete derivative/curl operators, and the r = 0 axis treatment | `src/grid.jl`, `src/operators.jl` | planned |
| [`pml.md`](pml.md) | complex coordinate stretching and the PML absorbing boundary | `src/operators.jl` | planned |
| [`materials.md`](materials.md) | ε/μ tensor construction, subpixel averaging, and Yee staggering | `src/materials.jl` | planned |
| [`system-assembly.md`](system-assembly.md) | how operators + materials form A x = b, sources, and BCs; eigenfrequency vs driven solves | `src/solver.jl` | drafted |
| [`postprocessing.md`](postprocessing.md) | Purcell factor, collection efficiency, far-field | `src/postprocess.jl` | planned |
| [`adjoint.md`](adjoint.md) | eigen-sensitivity gradients for optimization | `src/adjoint.jl` | planned |
| [`validation.md`](validation.md) | COMSOL and convergence validation methodology | `validation/` | planned |
| [`notation.md`](notation.md) | canonical symbol table | — | drafted |