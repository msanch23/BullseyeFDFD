# assets

Images referenced by the top-level `README.md`, and the two scripts that keep the README's
numbers and images in sync with the code. Keep images small (git stores every version forever).

## Updating the README

The README's report blocks and tables are generated from saved data; don't edit them by hand.

1. Re-create the data whose producer is affected by a change:
   - `julia --project=. assets/capture_example.jl` re-runs `examples/sanchezBullseye.jl` headless,
     saves its report text to `assets/readme/`, and saves `crosssection.png`, `eigenmode.png`,
     `farfield_kspace.png` and `driven_efield.png` here.
   - `validation/comsol/sweep_fdfd.jl`, then `validation/comsol/comsol_validation.ipynb`
     → `validation/comsol/results/summary.csv`.
   - `validation/adjoint/fd_vs_adjoint.ipynb` → `validation/adjoint/results/summary.csv`.
2. `julia assets/update_readme.jl` rewrites every block between
   `<!-- readme:PATH -->` and `<!-- /readme:PATH -->` from the file `PATH`
   (`.txt` → verbatim block, `.csv` → table). `--check` only reports whether any block is out of date.
   A table's headers and number formats are set by `HEADERS` and `FORMATS` in the script; a column
   with one value in every row (e.g. the grid) moves to a line under the table.

Every data file records the commit it came from. `update_readme.jl` warns when that commit's `src/`
or `Project.toml` differs from `HEAD`, i.e. when the numbers no longer describe the code.
