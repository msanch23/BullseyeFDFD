# assets

Static images referenced by the top-level `README.md`. Keep them small and
optimized (git stores every version forever).

## Regenerating `hero.png`

The hero image is a `show_sim` cross-section. To (re)generate it from the
reference design (see `examples/sanchezBullseye.jl`):

```julia
using BullseyeFDFD, Plots
# ... build design / substrates / geometry / grid / ϵ as in the example ...
savefig(show_sim(ϵ, grid), joinpath(@__DIR__, "hero.png"))
```
