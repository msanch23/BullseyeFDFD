using BullseyeFDFD
using Test
using SparseArrays
using LinearAlgebra

# ---------------------------------------------------------------------------
# Shared fixture: the sanchezBullseye reference case (matches examples/).
# Legacy hand-verified grid for these params was Nρ=1806, Nz=204; the current
# code gives 1805/202 — the small deltas are the intended stretch_grid
# sliver-merge + fill_layer→stretch_grid changes (stretched regions only).
# ---------------------------------------------------------------------------
const nClad = 1.00
const nSiN  = 2.01066
const nSiO2 = 1.45375
const nSi   = 3.69476 - 1im*0.00482
const tCBG  = 0.17886
const tSiO2 = 0.75057

const substrates = [(n = nSiO2, height = tSiO2),
                    (n = nSi,   height = Inf)]

const design = [0.355311, 0.100001, 0.382645, 0.100001, 0.362039, 0.100001,
                0.360574, 0.120593, 0.349761, 0.130700, 0.343483, 0.104798,
                0.367809, 0.100001, 0.379366, 0.102140, 0.377532, 0.114806,
                0.379070, 0.123844, 0.367272, 0.155886]

const λ      = 0.780
const gridkw = (; λPts = 32, geoPts = 32,
                  PMLρ = 0.5, PMLz_bot = 0.5, PMLz_top = 0.5,
                  padρ = 0.5, padz_bot = 0.5, padz_top = 0.5)

const geom = build_geometry(design, substrates, nClad, nSiN, tCBG)
const grid = conformal_grid(geom, λ; gridkw...)

@testset "BullseyeFDFD.jl" begin

    @testset "build_geometry" begin
        N_rings = (length(design) - 2) ÷ 2
        @test keys(geom) == (:eye, :trench, :ring, :buffer, :tCBG, :stack)
        @test length(geom.trench.width) == N_rings
        @test length(geom.ring.width)   == N_rings
        @test geom.eye.width    == design[1]        # eye = first design var
        @test geom.buffer.width == design[end]      # buffer = last design var
        @test geom.eye.n  == nSiN && geom.ring.n   == nSiN    # SiN ridges
        @test geom.trench.n == nClad && geom.buffer.n == nClad  # etched trenches
        @test geom.eye.height == tCBG && geom.trench.height == 0.0
        @test geom.tCBG == tCBG
        @test last(geom.stack) === geom.stack.wafer  # wafer is the last stack entry
        @test geom.stack.wafer.n == nSi

        # interleave mapping: pins trench[i]/ring[i] to the right design slots
        # (endpoint + length checks above wouldn't catch an off-by-one in the split)
        @test geom.trench.width == design[2:2:end-1]
        @test geom.ring.width   == design[3:2:end-1]

        # slotted path (idx_eye=2, prepends :slot). Prepending a slot to the same
        # design must leave eye/trench/ring/buffer identical to the unslotted case.
        dslot = [0.05; design]                      # length 23 ⇒ N=10
        gs = build_geometry(dslot, substrates, nClad, nSiN, tCBG; slotted=true)
        @test keys(gs) == (:slot, :eye, :trench, :ring, :buffer, :tCBG, :stack)
        @test gs.slot.n == nClad && gs.slot.height == 0.0 && gs.slot.width == 0.05
        @test gs.eye.width == design[1]             # eye shifts to design position 2
        @test gs.trench.width == geom.trench.width  # same rings/trenches as unslotted
        @test gs.ring.width   == geom.ring.width
        @test gs.buffer.width == geom.buffer.width

        # length guard (@assert) rejects a design that doesn't fit the topology
        @test_throws AssertionError build_geometry(design[1:end-1], substrates, nClad, nSiN, tCBG)
        @test_throws AssertionError build_geometry(dslot, substrates, nClad, nSiN, tCBG; slotted=false)

        # N substrates: dynamic substrate1..k naming + wafer
        subs3 = [(n = nSiO2, height = 0.30), (n = 1.60, height = 0.20), (n = nSi, height = Inf)]
        g3 = build_geometry(design, subs3, nClad, nSiN, tCBG)
        @test keys(g3.stack) == (:substrate1, :substrate2, :wafer)
        @test g3.stack.substrate1.n == nSiO2 && g3.stack.substrate1.height == 0.30
        @test g3.stack.substrate2.n == 1.60
        @test g3.stack.wafer.height == Inf

        # wafer-only stack (k=0 edge case: empty substrate NamedTuple)
        g1 = build_geometry(design, [(n = nSi, height = Inf)], nClad, nSiN, tCBG)
        @test keys(g1.stack) == (:wafer,)
    end

    @testset "conformal_grid — array-length invariants" begin
        @test length(grid.ρp) == grid.Nρ + 1        # node arrays are N+1
        @test length(grid.zp) == grid.Nz + 1
        @test length(grid.Δρp) == grid.Nρ           # width / dual arrays are N
        @test length(grid.Δzp) == grid.Nz
        @test length(grid.ρd) == grid.Nρ
        @test length(grid.zd) == grid.Nz
        @test length(grid.Δρd) == grid.Nρ
        @test length(grid.Δzd) == grid.Nz
        @test grid.N == grid.Nρ * grid.Nz
    end

    @testset "conformal_grid — coordinate sanity" begin
        @test grid.ρp[1] == 0.0
        @test grid.zp[1] == 0.0
        @test issorted(grid.ρp)
        @test issorted(grid.zp)
        @test all(grid.Δρp .> 0) && all(grid.Δzp .> 0)
        @test all(grid.Δρd .> 0) && all(grid.Δzd .> 0)
        @test grid.R ≈ grid.ρp[end]
        @test grid.Z ≈ grid.zp[end]
    end

    @testset "conformal_grid — physical/index consistency" begin
        # device layer spans exactly tCBG and lands on a cell edge
        @test grid.cbg_top - grid.cbg_bot ≈ tCBG
        @test grid.zp[grid.idx_cbg_bot] ≈ grid.cbg_bot
        # domain = features + pad + PML on each side
        @test grid.R - grid.total_cbg_radius ≈ gridkw.padρ + gridkw.PMLρ
        @test grid.Z - grid.cbg_top ≈ gridkw.padz_top + gridkw.PMLz_top
        # index bookkeeping is in-range and ordered
        @test 1 ≤ grid.idx_cbg_bot ≤ grid.idx_cbg_cen ≤ grid.idx_cbg_top ≤ grid.Nz
        @test 1 ≤ grid.idx_PMLρ ≤ grid.Nρ
        @test 1 ≤ grid.idx_PMLz_bot ≤ grid.idx_PMLz_top ≤ grid.Nz
    end

    # Exact cell counts — a characterization check. These SHOULD change if you
    # intentionally retune the meshing (resolution, stretch ratios, pad/PML);
    # update the numbers when you do. A silent change here means a bug.
    @testset "conformal_grid — characterization (cell counts)" begin
        @test grid.Nρ == 1805
        @test grid.Nz == 202
        @test grid.N  == 364610
    end

    # -- materials: build_epsilon ------------------------------------------
    ε_ref = build_epsilon(geom, grid)

    @testset "build_epsilon — shape & smoke" begin
        @test ε_ref isa Matrix{ComplexF64}
        @test size(ε_ref) == (grid.Nρ, grid.Nz)
        @test !any(isnan, ε_ref)
        @test all(isfinite, ε_ref)
    end

    @testset "build_epsilon — material values" begin
        jc = grid.idx_cbg_cen
        @test ε_ref[1, jc]   ≈ ComplexF64(nSiN)^2               # eye disk = SiN
        @test ε_ref[end, jc] ≈ ComplexF64(nSiN)^2               # bulk SiN beyond the buffer
        band = real.(ε_ref[:, jc])
        @test any(x -> isapprox(x, nClad^2; atol=1e-9), band)   # etched trenches = clad
        @test any(x -> isapprox(x, nSiN^2;  atol=1e-6), band)   # SiN ridges
        @test ε_ref[1, grid.idx_cbg_bot - 1] ≈ ComplexF64(nSiO2)^2  # top substrate
        @test ε_ref[1, 1]   ≈ ComplexF64(nSi)^2                     # wafer (semi-infinite)
        @test imag(ε_ref[1, 1]) < 0                                 # absorptive
        @test ε_ref[1, end] ≈ ComplexF64(nClad)^2                   # cladding above
    end

    # Slotted cavity: a central clad hole (ρ=0) precedes the SiN eye. build_epsilon
    # must prepend the slot to match grid.jl's slot→eye→… radial layout.
    @testset "build_epsilon — slotted (central clad hole)" begin
        dslot = [0.05; design]                          # prepend slot ⇒ N=10, length 23
        gs  = build_geometry(dslot, substrates, nClad, nSiN, tCBG; slotted = true)
        grs = conformal_grid(gs, λ; gridkw...)
        εs  = build_epsilon(gs, grs)
        jc  = grs.idx_cbg_cen
        @test εs[1, jc] ≈ ComplexF64(nClad)^2           # slot = clad hole at ρ=0, NOT the SiN eye
        ie = findfirst(i -> grs.ρd[i] > gs.slot.width + 0.5*gs.eye.width, 1:grs.Nρ)
        @test εs[ie, jc] ≈ ComplexF64(nSiN)^2           # eye ring just outside the slot
        @test grs.total_cbg_radius ≈ gs.slot.width + gs.eye.width +
              sum(gs.trench.width) + sum(gs.ring.width) + gs.buffer.width   # ε aligns with grid
    end

    # All indices equal ⇒ ε uniform everywhere: proves the ρ- and z-subpixel
    # fractions form a partition of unity (no material dropped/double-counted).
    @testset "build_epsilon — partition (uniform index ⇒ uniform ε)" begin
        n0 = 1.7
        subs0 = [(n = n0, height = tSiO2), (n = n0, height = Inf)]
        geom0 = build_geometry(design, subs0, n0, n0, tCBG)
        grid0 = conformal_grid(geom0, λ; gridkw...)
        ε0 = build_epsilon(geom0, grid0)
        @test all(x -> isapprox(x, ComplexF64(n0)^2; atol=1e-12), ε0)
    end

    # Guards the adjoint-critical property: on a FROZEN grid, ε varies smoothly
    # with widths (a mid-cell interface yields a blended cell) and thicknesses.
    @testset "build_epsilon — subpixel active on a frozen grid" begin
        jc = grid.idx_cbg_cen
        gd = build_geometry([design[1] + 0.007; design[2:end]], substrates, nClad, nSiN, tCBG)
        εd = build_epsilon(gd, grid)                            # SAME (frozen) grid
        @test any(x -> nClad^2 + 1e-6 < real(x) < nSiN^2 - 1e-6, εd[:, jc])  # blended ρ-cell
        Ndev = grid.idx_cbg_top - grid.idx_cbg_bot + 1
        gt = build_geometry(design, substrates, nClad, nSiN, tCBG + 0.5*tCBG/Ndev)
        @test build_epsilon(gt, grid) != ε_ref                 # ε responds to tCBG on a frozen grid
    end

    # -- operators: CeChOperators ------------------------------------------
    # A full-size smoke test, then structural/analytic checks on small grids
    # where slicing the 3N×3N blocks is cheap. The block layout is
    #   [ (1,1) (1,2) (1,3) ;  (2,1) (2,2) (2,3) ;  (3,1) (3,2) (3,3) ]
    # each block N×N; (1,2) = -Dz_E_2D, (1,3) & (3,1) carry the im·m/ρ metric.
    @testset "CeChOperators — smoke (full grid)" begin
        Ce, Ch = CeChOperators(grid, λ)
        n3 = 3 * grid.N
        @test size(Ce) == (n3, n3) && size(Ch) == (n3, n3)
        @test eltype(Ce) == ComplexF64 && eltype(Ch) == ComplexF64
        @test all(isfinite, nonzeros(Ce)) && !any(isnan, nonzeros(Ce))
        @test all(isfinite, nonzeros(Ch)) && !any(isnan, nonzeros(Ch))
    end

    # coarse grids (fast to slice); one with PML, one without
    gpml = conformal_grid(geom, λ; λPts = 6, geoPts = 2,
                          PMLρ = 0.2, PMLz_bot = 0.2, PMLz_top = 0.2,
                          padρ = 0.2, padz_bot = 0.2, padz_top = 0.2)
    g0   = conformal_grid(geom, λ; λPts = 6, geoPts = 2,
                          PMLρ = 0.0, PMLz_bot = 0.0, PMLz_top = 0.0,
                          padρ = 0.2, padz_bot = 0.2, padz_top = 0.2)

    @testset "CeChOperators — block structure" begin
        N = gpml.N
        Ce, Ch = CeChOperators(gpml, λ; m = 1)
        for B in (Ce, Ch), b in 0:2                 # the on-diagonal blocks are zero
            @test iszero(B[b*N+1:(b+1)*N, b*N+1:(b+1)*N])
        end
    end

    @testset "CeChOperators — azimuthal (m) coupling" begin
        N = gpml.N
        Ce0, Ch0 = CeChOperators(gpml, λ; m = 0)    # ρ↔z metric blocks vanish at m=0 ...
        @test iszero(Ce0[1:N, 2N+1:3N]) && iszero(Ce0[2N+1:3N, 1:N])
        @test iszero(Ch0[1:N, 2N+1:3N]) && iszero(Ch0[2N+1:3N, 1:N])
        Ce1, _ = CeChOperators(gpml, λ; m = 1)      # ... and return for m ≠ 0
        @test !iszero(Ce1[1:N, 2N+1:3N])
    end

    @testset "CeChOperators — PML stretch path" begin
        # m = 0 removes the imaginary metric terms, so any imaginary part is PML
        Ce_p, _ = CeChOperators(gpml, λ; m = 0)
        @test any(!iszero, imag.(nonzeros(Ce_p)))   # stretch is active with PML on
        Ce0, Ch0 = CeChOperators(g0, λ; m = 0)
        @test all(iszero, imag.(nonzeros(Ce0)))     # PML off ⇒ purely real
        @test all(iszero, imag.(nonzeros(Ch0)))
    end

    @testset "CeChOperators — z-derivative stencil" begin
        Ce, _ = CeChOperators(g0, λ; m = 0)         # PML off ⇒ real derivative stencil
        N, Nρ, Nz = g0.N, g0.Nρ, g0.Nz
        Dz = Ce[1:N, N+1:2N]                        # block (1,2) = -Dz_E_2D

        slope = 0.37                                # a field linear in z: f = slope·z
        f = [slope * g0.zp[j] for j in 1:Nz for i in 1:Nρ]   # ρ fastest, z slowest
        r = Dz * f                                  # ∂z of a ramp ⇒ -slope, interior only

        interior = [(j-1)*Nρ + i for j in 1:Nz-1 for i in 1:Nρ]  # drop last (one-sided) z-row
        @test all(x -> isapprox(x, -slope; atol = 1e-10), real.(r[interior]))
        @test all(iszero, imag.(r[interior]))
    end

    @testset "CeChOperators — ρ-derivative stencil" begin
        Ce, _ = CeChOperators(g0, λ; m = 0)         # PML off ⇒ real derivative stencil
        N, Nρ, Nz = g0.N, g0.Nρ, g0.Nz
        Dρ = Ce[N+1:2N, 2N+1:3N]                    # block (2,3) = -Dρ_E_2D

        slope = 0.29                                # a field linear in ρ: f = slope·ρ
        f = [slope * g0.ρp[i] for j in 1:Nz for i in 1:Nρ]   # ρ fastest, z slowest
        r = Dρ * f                                  # ∂ρ of a ramp ⇒ -slope, interior only

        interior = [(j-1)*Nρ + i for j in 1:Nz for i in 1:Nρ-1]  # drop last (one-sided) ρ-row
        @test all(x -> isapprox(x, -slope; atol = 1e-10), real.(r[interior]))
        @test all(iszero, imag.(r[interior]))
    end

    # The ρ=0 axis is the cylindrical-specific singularity: 1/ρ is forced to 0
    # there (the im·m/ρ metric block must vanish on-axis, else the curl blows up).
    @testset "CeChOperators — ρ=0 axis handling" begin
        Ce, _ = CeChOperators(gpml, λ; m = 1)
        N, Nρ, Nz = gpml.N, gpml.Nρ, gpml.Nz
        d = diag(Ce[1:N, 2N+1:3N])                  # (1,3) block = -im·m·(1/ρp), diagonal
        axis = [(j-1)*Nρ + 1 for j in 1:Nz]         # ρ = 0 nodes (i = 1)
        @test all(iszero, d[axis])                  # 1/ρp = 0 on the axis
        @test all(!iszero, d[setdiff(1:N, axis)])   # nonzero everywhere off-axis
    end

    # -- materials: Yee-staggered inverse tensors --------------------------
    @testset "material_tensors" begin
        invϵ, invμ = BullseyeFDFD.material_tensors(ε_ref, grid)
        N = grid.N
        @test size(invϵ) == (3N, 3N) && size(invμ) == (3N, 3N)
        @test eltype(invϵ) == ComplexF64
        @test nnz(invϵ) == 3N && nnz(invμ) == 3N            # strictly diagonal
        @test all(≈(1.0 + 0im), diag(invμ))                 # μ = I (non-magnetic)
        @test all(isfinite, diag(invϵ))
        # Ez block (3rd) diag = 1/ϵ_zz with ϵ_zz[i,j] = ½(ϵ[i,j] + ϵ[i-1,j])
        i = grid.Nρ ÷ 2; j = grid.idx_cbg_cen
        @test diag(invϵ)[2N + (j-1)*grid.Nρ + i] ≈ 1 / (0.5 * (ε_ref[i,j] + ε_ref[i-1,j]))
        # uniform index ⇒ every staggered tensor = n0² ⇒ diag = 1/n0²
        n0 = 1.7; subs0 = [(n = n0, height = tSiO2), (n = n0, height = Inf)]
        g0u  = build_geometry(design, subs0, n0, n0, tCBG)
        gr0u = conformal_grid(g0u, λ; gridkw...)
        iϵ0, _ = BullseyeFDFD.material_tensors(build_epsilon(g0u, gr0u), gr0u)
        @test all(x -> isapprox(x, 1 / ComplexF64(n0)^2; atol = 1e-10), diag(iϵ0))
    end

    # -- solver: solve_sim + dipole ----------------------------------------
    # Dedicated coarse grid — eig (Arpack), driven (LU), Purcell bulk solve, and
    # the far-field crunch must all run cheaply. report/plots off = quiet & headless.
    gsolve = conformal_grid(geom, λ; λPts=10, geoPts=6,
                            PMLρ=0.3, PMLz_bot=0.3, PMLz_top=0.3,
                            padρ=0.3, padz_bot=0.3, padz_top=0.3)
    εsolve = build_epsilon(geom, gsolve)

    @testset "dipole — constructor (grid-first)" begin
        d0 = dipole(gsolve)
        @test d0.λ === nothing && d0.pol === :ρ && d0.ρ_idx == 2 && d0.z_idx == gsolve.idx_cbg_cen
        @test dipole(gsolve, 0.781).λ == 0.781           # grid-first positional λ
        @test dipole(gsolve; pol=:z, ρ_idx=5).pol === :z
    end

    @testset "solve_sim — source selection rule (pol↔m)" begin
        @test BullseyeFDFD._check_source((; pol=:z), 0)      # ẑ ⇒ m=0
        @test BullseyeFDFD._check_source((; pol=:ρ), 1)      # in-plane ⇒ m=1
        @test_throws ErrorException BullseyeFDFD._check_source((; pol=:z), 1)
    end

    @testset "solve_sim — eigenmodes" begin
        sim = solve_sim(gsolve, εsolve, λ; Nmodes=2, report=false, plots=false)
        @test length(sim.modes) == 2
        @test all(m -> m.λ > 0 && isfinite(m.λ) && !isnan(m.Q), sim.modes)
        @test all(m -> all(isfinite, m.Eρ), sim.modes)
        @test all(m -> hasproperty(m, :CF) && hasproperty(m, :Fp), sim.modes)
        @test sim.driven === nothing
        @test size(sim.Ce) == (3gsolve.N, 3gsolve.N)         # Ce surfaced for postprocess
    end

    @testset "solve_sim — driven dipole" begin
        sd = solve_sim(gsolve, εsolve, λ; source=dipole(gsolve, 0.781),
                       report=false, plots=false, purcell=false)
        @test sd.modes === nothing
        @test sd.driven.λ == 0.781                           # explicit drive λ honored
        @test all(isfinite, sd.driven.raw_E)
        @test isnan(sd.driven.Purcell)                       # purcell=false ⇒ not computed
        se = solve_sim(gsolve, εsolve, λ; Nmodes=2, source=dipole(gsolve),
                       report=false, plots=false, purcell=false)
        @test se.driven.λ == argmax(m -> m.CF, se.modes).λ   # dipole(grid) ⇒ CF-selected resonance
        sp = solve_sim(gsolve, εsolve, λ; Nmodes=1, source=dipole(gsolve),
                       purcell=true, report=false, plots=false)
        @test isfinite(sp.driven.Purcell) && sp.driven.Purcell > 0   # LDOS bulk ratio
    end

    @testset "solve_sim — NA gates collection FOMs" begin
        s0 = solve_sim(gsolve, εsolve, λ;      Nmodes=1, report=false, plots=false)
        @test !hasproperty(s0.modes[1], :η)                  # NA omitted ⇒ no far-field crunch
        s1 = solve_sim(gsolve, εsolve, λ, 0.4; Nmodes=1, report=false, plots=false)
        @test hasproperty(s1.modes[1], :η) && isfinite(s1.modes[1].η)
    end
end
