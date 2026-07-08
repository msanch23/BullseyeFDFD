# Theory notes

Derivations and background for the cylindrical FDFD solver, versioned alongside
the code they describe. One file per topic; keep each note close to the `src/`
file it explains.

Suggested topics:

- `maxwell-cylindrical.md` — Maxwell's equations in (r, φ, z), azimuthal mode
  decomposition (e^{imφ}), and the reduction to a 2D (r, z) problem.
- `discretization.md` — Yee-cell layout in cylindrical coordinates, the discrete
  derivative/curl operators, and the r = 0 axis treatment.
- `pml.md` — complex coordinate stretching and the PML absorbing boundary.
- `system-assembly.md` — how operators + materials form A x = b, sources, and BCs.
- `postprocessing.md` — Purcell factor, collection efficiency, far-field.

LaTeX in `$…$` / `$$…$$` renders on GitHub. Cross-link each note to its
`src/` counterpart and to the `examples/` script that demonstrates it.
