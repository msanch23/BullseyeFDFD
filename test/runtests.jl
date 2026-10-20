using BullseyeFDFD
using Test
using SparseArrays
using LinearAlgebra

# ---------------------------------------------------------------------------
# Shared fixture: the sanchezBullseye reference case (matches examples/).
# Legacy hand-verified grid for these params was Nρ=1806, Nz=204; the
# stretch_grid sliver-merge + fill_layer→stretch_grid changes gave 1805/202,
# the material-aware substrate/wafer z mesh (fit_layer) gave 1805/223, and tying
# the pad cap to λPts (1.5λ/(n·λPts)) gives 1805/248.
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
        @test grid.Nz == 248
        @test grid.N  == 447640
    end

    # The substrate grades down to the wafer's λ/(n·λPts) at the interface and the wafer starts
    # there; before, stretch_grid's remainder cell sat on the interface and set the first Si cell.
    @testset "conformal_grid — substrate → wafer interface" begin
        z_if = grid.cbg_bot - tSiO2
        i    = argmin(abs.(grid.zp .- z_if))
        Δw   = λ / (real(nSi) * gridkw.λPts)
        @test grid.zp[i] ≈ z_if                     # interface on a cell edge
        @test grid.Δzp[i]   ≈ Δw                    # last SiO2 cell
        @test grid.Δzp[i-1] ≈ Δw                    # first Si cell
        @test all(grid.Δzp[2:i-1] .≤ 1.5λ / (real(nSi) * gridkw.λPts) + 1e-12)   # wafer cap 1.5λ/(n·λPts)
        for dt in (-2e-4, 2e-4)                     # no remainder flip under a small tSiO2 change
            gt = conformal_grid(build_geometry(design, [(n = nSiO2, height = tSiO2 + dt), (n = nSi, height = Inf)],
                                               nClad, nSiN, tCBG), λ; gridkw...)
            @test gt.Nz == grid.Nz
        end
    end

    # Top pad mirrors the bottom: capped by the cladding's 1.5λ/(n·λPts). A SiO2 cladding at
    # geoPts 8 makes that cap bind (3·target_Δz is coarser there).
    @testset "conformal_grid — material-aware top pad" begin
        gc  = build_geometry(design, substrates, nSiO2, nSiN, tCBG)
        grc = conformal_grid(gc, λ; gridkw..., λPts = 24, geoPts = 8)
        top = grc.Δzp[grc.idx_cbg_top+1:end-1]      # the last cell may absorb stretch_grid's remainder
        cap = 1.5λ / (nSiO2 * 24)
        @test maximum(top) ≤ cap + 1e-12
        @test maximum(top) > 0.99 * cap
    end

    @testset "fit_layer" begin
        c = BullseyeFDFD.fit_layer(0.75, 0.0056, 0.0066, 1.15, 0.018)
        @test sum(c) ≈ 0.75
        @test c[1] ≈ 0.0056 && c[end] ≈ 0.0066
        @test all(c .≤ 0.018)
        t = BullseyeFDFD.fit_layer(0.03, 0.0112, nothing, 1.15, 0.022)     # too thin to grade
        @test sum(t) ≈ 0.03 && all(t .≈ t[1])
    end

    # -- materials: build_epsilon ------------------------------------------
    ε_ref = build_epsilon(geom, grid)

    @testset "build_epsilon — shape & smoke" begin
        @test ε_ref isa NamedTuple{(:cell, :ρρ, :ϕϕ, :zz)}
        for M in ε_ref
            @test M isa Matrix{ComplexF64}
            @test size(M) == (grid.Nρ, grid.Nz)
            @test all(isfinite, M)
        end
    end

    @testset "build_epsilon — material values" begin
        jc = grid.idx_cbg_cen
        @test ε_ref.cell[1, jc]   ≈ ComplexF64(nSiN)^2               # eye disk = SiN
        @test ε_ref.cell[end, jc] ≈ ComplexF64(nSiN)^2               # bulk SiN beyond the buffer
        band = real.(ε_ref.cell[:, jc])
        @test any(x -> isapprox(x, nClad^2; atol=1e-9), band)   # etched trenches = clad
        @test any(x -> isapprox(x, nSiN^2;  atol=1e-6), band)   # SiN ridges
        @test ε_ref.cell[1, grid.idx_cbg_bot - 1] ≈ ComplexF64(nSiO2)^2  # top substrate
        @test ε_ref.cell[1, 1]   ≈ ComplexF64(nSi)^2                     # wafer (semi-infinite)
        @test imag(ε_ref.cell[1, 1]) < 0                                 # absorptive
        @test ε_ref.cell[1, end] ≈ ComplexF64(nClad)^2                   # cladding above
        # away from every wall, all three Yee components see the cell's material
        @test ε_ref.ρρ[1, jc] ≈ ε_ref.ϕϕ[1, jc] ≈ ε_ref.zz[1, jc] ≈ ε_ref.cell[1, jc]
    end

    # Slotted cavity: a central clad hole (ρ=0) precedes the SiN eye. build_epsilon
    # must prepend the slot to match grid.jl's slot→eye→… radial layout.
    @testset "build_epsilon — slotted (central clad hole)" begin
        dslot = [0.05; design]                          # prepend slot ⇒ N=10, length 23
        gs  = build_geometry(dslot, substrates, nClad, nSiN, tCBG; slotted = true)
        grs = conformal_grid(gs, λ; gridkw...)
        εs  = build_epsilon(gs, grs)
        jc  = grs.idx_cbg_cen
        @test εs.cell[1, jc] ≈ ComplexF64(nClad)^2      # slot = clad hole at ρ=0, NOT the SiN eye
        @test εs.zz[1, jc]   ≈ ComplexF64(nClad)^2      # Ez on the axis sits in the slot too
        ie = findfirst(i -> grs.ρd[i] > gs.slot.width + 0.5*gs.eye.width, 1:grs.Nρ)
        @test εs.cell[ie, jc] ≈ ComplexF64(nSiN)^2      # eye ring just outside the slot
        @test grs.total_cbg_radius ≈ gs.slot.width + gs.eye.width +
              sum(gs.trench.width) + sum(gs.ring.width) + gs.buffer.width   # ε aligns with grid
    end

    # All indices equal ⇒ ε uniform everywhere: proves the ρ- and z-subpixel
    # fractions (ramped walls included) form a partition of unity in every map.
    @testset "build_epsilon — partition (uniform index ⇒ uniform ε)" begin
        n0 = 1.7
        subs0 = [(n = n0, height = tSiO2), (n = n0, height = Inf)]
        geom0 = build_geometry(design, subs0, n0, n0, tCBG)
        grid0 = conformal_grid(geom0, λ; gridkw...)
        ε0 = build_epsilon(geom0, grid0)
        @test all(M -> all(x -> isapprox(x, ComplexF64(n0)^2; atol=1e-12), M), ε0)
    end

    # Guards the adjoint-critical property: on a FROZEN grid, ε varies smoothly
    # with widths (a mid-cell interface yields a blended cell) and thicknesses.
    @testset "build_epsilon — subpixel active on a frozen grid" begin
        jc = grid.idx_cbg_cen
        gd = build_geometry([design[1] + 0.007; design[2:end]], substrates, nClad, nSiN, tCBG)
        εd = build_epsilon(gd, grid)                            # SAME (frozen) grid
        @test any(x -> nClad^2 + 1e-6 < real(x) < nSiN^2 - 1e-6, εd.cell[:, jc])  # blended ρ-cell
        Ndev = grid.idx_cbg_top - grid.idx_cbg_bot + 1
        gt = build_geometry(design, substrates, nClad, nSiN, tCBG + 0.5*tCBG/Ndev)
        @test build_epsilon(gt, grid) != ε_ref                 # ε responds to tCBG on a frozen grid
    end

    # With sharp walls (ramp = (0, 0)), a wall inside a voxel gives the harmonic mean for the
    # field component normal to it and the area/length-weighted arithmetic mean for the others.
    @testset "build_epsilon — anisotropic averaging (sharp walls)" begin
        eS, eC = ComplexF64(nSiN)^2, ComplexF64(nClad)^2
        jc = grid.idx_cbg_cen
        i0 = argmin(abs.(grid.ρp .- design[1]))                 # eye wall's node at the seed
        ρw = design[1] + 0.3 * grid.Δρp[i0]                     # 0.3 cell outward: inside cell i0
        εd = build_epsilon(build_geometry([ρw; design[2:end]], substrates, nClad, nSiN, tCBG), grid;
                           ramp = (0.0, 0.0))
        a, b = grid.ρp[i0], grid.ρp[i0+1]; f = (ρw^2 - a^2) / (b^2 - a^2)
        @test εd.ρρ[i0, jc]   ≈ 1 / (f/eS + (1 - f)/eC)          # Eρ ⟂ ρ-wall: harmonic
        @test εd.cell[i0, jc] ≈ f*eS + (1 - f)*eC                # cell map: arithmetic
        a, b = grid.ρd[i0-1], grid.ρd[i0]; g = (ρw^2 - a^2) / (b^2 - a^2)
        @test εd.ϕϕ[i0, jc] ≈ g*eS + (1 - g)*eC                  # Eϕ, Ez ∥ ρ-wall: arithmetic over
        @test εd.zz[i0, jc] ≈ g*eS + (1 - g)*eC                  # the node voxel [ρd[i0-1], ρd[i0]]
        j0 = grid.idx_cbg_top + 1                               # first cladding cell at the seed
        zw = grid.cbg_top + 0.3 * grid.Δzp[j0]                  # top wall 0.3 cell up: inside cell j0
        εt = build_epsilon(build_geometry(design, substrates, nClad, nSiN, zw - grid.cbg_bot), grid;
                           ramp = (0.0, 0.0))
        @test εt.zz[1, j0]   ≈ 1 / (0.3/eS + 0.7/eC)             # Ez ⟂ z-wall: harmonic
        @test εt.cell[1, j0] ≈ 0.3*eS + 0.7*eC
        h = (zw - grid.zd[j0-1]) / (grid.zd[j0] - grid.zd[j0-1])
        @test εt.ρρ[1, j0] ≈ h*eS + (1 - h)*eC                   # Eρ, Eϕ ∥ z-wall: arithmetic over
        @test εt.ϕϕ[1, j0] ≈ h*eS + (1 - h)*eC                   # the node voxel [zd[j0-1], zd[j0]]
    end

    # The ramped walls make the Yee-node ε C¹ in every width and thickness. At the conformal seed
    # every wall sits on a grid line, where the sharp map kinks; one-sided derivatives must agree.
    @testset "build_epsilon — C¹ on a frozen grid" begin
        yee(d, t, s; kw...) = BullseyeFDFD._yee_eps(build_epsilon(build_geometry(d,
                                  [(n = nSiO2, height = s), (n = nSi, height = Inf)], nClad, nSiN, t), grid; kw...), grid)
        nd = length(design)
        bump(k, h) = k <= nd ? (setindex!(copy(design), design[k] + h, k), tCBG, tSiO2) :
                     k == nd + 1 ? (design, tCBG + h, tSiO2) : (design, tCBG, tSiO2 + h)
        δ = 1e-7; y0 = yee(design, tCBG, tSiO2)
        for k in (1, 2, 3, nd + 1, nd + 2)                      # eye, trench 1, ring 1, tCBG, tSiO2
            R = (yee(bump(k, δ)...) .- y0) ./ δ; L = (y0 .- yee(bump(k, -δ)...)) ./ δ
            @test maximum(abs.(R .- L)) < 1e-3 * maximum(abs.(R))
        end
        s0 = yee(design, tCBG, tSiO2; ramp = (0.0, 0.0))          # control: sharp walls DO kink here
        R = (yee(bump(1, δ)...; ramp = (0.0, 0.0)) .- s0) ./ δ
        L = (s0 .- yee(bump(1, -δ)...; ramp = (0.0, 0.0))) ./ δ
        @test maximum(abs.(R .- L)) > 0.1 * maximum(abs.(R))
    end

    # Default ramp widths come from the frozen grid. On the default mesh every z wall has equal
    # cells on both sides, and each wall gets exactly its own cell (a hat kernel). The slab's walls
    # no longer take the substrate/wafer cell, which made them depend on the wafer index.
    @testset "build_epsilon — default ramp widths" begin
        gd = conformal_grid(geom, λ)                            # default 24/16 mesh
        r  = BullseyeFDFD._default_ramp(gd)
        jb, jt = gd.idx_cbg_bot, gd.idx_cbg_top
        js = argmin(abs.(gd.zp .- (gd.cbg_bot - tSiO2)))        # node of the SiO2/Si wall
        @test r.slab ≈ gd.Δzp[jt] ≈ gd.Δzp[jt+1] ≈ gd.Δzp[jb-1] # slab walls: their own cell
        @test r.substrate ≈ gd.Δzp[js-1] ≈ gd.Δzp[js]           # SiO2/Si wall: its own cell
        @test r.substrate < r.slab                              # one shared width gave the slab this
        gq = conformal_grid(build_geometry(design, [(n = nSiO2, height = tSiO2), (n = nSiO2, height = Inf)],
                                           nClad, nSiN, tCBG), λ)  # same stack on a SiO2 wafer
        @test BullseyeFDFD._default_ramp(gq).slab == r.slab      # slab walls ignore the wafer index
        # the legacy pair (wρ, wz) means one width for every z wall
        e2 = build_epsilon(geom, grid; ramp = (r.ρ, r.slab))
        e3 = build_epsilon(geom, grid; ramp = (ρ = r.ρ, slab = r.slab, substrate = r.slab))
        @test all(getfield(e2, f) == getfield(e3, f) for f in keys(e2))
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
        N, Nρ, Nz = gpml.N, gpml.Nρ, gpml.Nz
        Ce, Ch = CeChOperators(gpml, λ; m = 1)
        ax = [(j-1)*Nρ + 1 for j in 1:Nz]           # ρ = 0 rows within a block
        for B in (Ce, Ch), b in 0:2                 # the on-diagonal blocks are zero ...
            D = B[b*N+1:(b+1)*N, b*N+1:(b+1)*N]
            B === Ch && b == 1 && (D = D[setdiff(1:N, ax), :])
            @test iszero(D)
        end
        @test !iszero(Ch[N .+ ax, N+1:2N])          # ... except Ch's axis Eϕ rows: Hρ(0) = i·m·Hϕ[1]
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

    # Radial PML: ρ → ρ̃ = ∫₀^ρ S dρ′ stretches the derivatives AND every explicit ρ.
    @testset "s_vectors — stretched coordinate ρ̃" begin
        p, R_target = 3.5, 1e-12
        _, _, _, _, ρ̃p, ρ̃d = BullseyeFDFD.s_vectors(gpml, λ; p, R_target)
        ρ_start = gpml.R - gpml.PMLρ
        ρp = gpml.ρp[1:end-1]
        @test ρ̃p[ρp .<= ρ_start] == ρp[ρp .<= ρ_start]               # ρ̃ ≡ ρ outside the layer
        @test ρ̃d[gpml.ρd .<= ρ_start] == gpml.ρd[gpml.ρd .<= ρ_start]
        L, u = gpml.PMLρ, gpml.ρd[end] - ρ_start                      # closed form at the last node
        σ_max = -(p + 1) * log(R_target) / (2L)
        @test ρ̃d[end] ≈ gpml.ρd[end] - im * (σ_max * λ / 2π) * L / (p + 1) * (u / L)^(p + 1)
        _, _, _, _, ρ̃p0, ρ̃d0 = BullseyeFDFD.s_vectors(g0, λ)
        @test ρ̃p0 == g0.ρp[1:end-1] && ρ̃d0 == g0.ρd                  # PMLρ = 0 ⇒ ρ̃ = ρ

        # dρ̃/dρ = S: ρ̃ and the derivative stretch describe the same coordinate
        _, Sρd, _, _, ρ̃pf, _ = BullseyeFDFD.s_vectors(grid, λ)
        i = findall(grid.ρp[2:end-1] .> grid.R - grid.PMLρ)
        @test all(isapprox.(diff(ρ̃pf)[i] ./ grid.Δρp[i], Sρd[i]; rtol = 1e-2))
    end

    @testset "CeChOperators — PML metric uses ρ̃" begin
        Ce, Ch = CeChOperators(gpml, λ; m = 1)
        N, Nρ, Nz = gpml.N, gpml.Nρ, gpml.Nz
        _, _, _, _, ρ̃p, ρ̃d = BullseyeFDFD.s_vectors(gpml, λ)
        inpml = repeat(gpml.ρp[1:end-1] .> gpml.R - gpml.PMLρ, Nz)
        d = diag(Ce[1:N, 2N+1:3N])                  # (1,3) block = -im·m/ρ̃p: real part only in the PML
        @test all(!iszero, real.(d[inpml])) && all(iszero, real.(d[.!inpml]))

        # ρ̃Eϕ = const ⇒ ρ̃⁻¹∂ρ(ρ̃Eϕ) = 0 on interior rows, PML included (inner AND outer ρ̃)
        e = [i == 1 ? 0.0im : 1 / ρ̃p[i] for j in 1:Nz for i in 1:Nρ]
        r = Ce[2N+1:3N, N+1:2N] * e                 # (3,2) block = Dρ_ρE
        @test all(x -> abs(x) < 1e-8, r[[(j-1)*Nρ + i for j in 1:Nz for i in 2:Nρ-1]])
        h = [1 / ρ̃d[i] for j in 1:Nz for i in 1:Nρ]
        r = Ch[2N+1:3N, N+1:2N] * h                 # (3,2) block = Dρ_ρH
        @test all(x -> abs(x) < 1e-8, r[[(j-1)*Nρ + i for j in 1:Nz for i in 2:Nρ]])
    end

    # ρ = 0 at |m| = 1: the axis Eϕ row is −i·m × the Eρ[1] row (Hρ(0) = i·m·Hϕ(0)) and the
    # axis Hρ keeps −i·m·Ez[2]/ρp[2] (Seet A18), so no free 1D chain lives on the axis.
    @testset "CeChOperators — ρ=0 axis regularity rows" begin
        N, Nρ, Nz = gpml.N, gpml.Nρ, gpml.Nz
        r0 = [(j-1)*Nρ + 1 for j in 1:Nz]           # Eρ[1] / Hρ(0) rows
        rϕ = N .+ r0                                # Eϕ(0) rows
        for m in (1, -1)
            Ce, Ch = CeChOperators(gpml, λ; m)
            @test Ch[rϕ, :] == -1im * m * Ch[r0, :]
            @test all(j -> Ce[r0[j], 2N + r0[j] + 1] ≈ -1im * m / gpml.ρp[2], 1:Nz)
            A = Ch * Ce                             # ε = μ = 1
            @test iszero(A[setdiff(1:3N, rϕ), rϕ])  # nothing reads the axis Eϕ ...
            @test iszero(A[rϕ, rϕ])                 # ... and it has no dynamics of its own
        end
        for m in (0, 2)
            _, Ch = CeChOperators(gpml, λ; m)
            @test iszero(Ch[rϕ, :])                 # Eϕ(0) = 0 for |m| ≠ 1
        end
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
        # build_epsilon output: diag = 1 ./ [ρρ; ϕϕ; zz], each block column-major (ρ fastest)
        @test diag(invϵ) == 1 ./ vcat(vec(ε_ref.ρρ), vec(ε_ref.ϕϕ), vec(ε_ref.zz))
        # plain cell map (hand-built or uniform): ½/½ Yee averaging, e.g. Ez: ½(ϵ[i,j] + ϵ[i-1,j])
        invc, _ = BullseyeFDFD.material_tensors(ε_ref.cell, grid)
        i = grid.Nρ ÷ 2; j = grid.idx_cbg_cen
        @test diag(invc)[2N + (j-1)*grid.Nρ + i] ≈ 1 / (0.5 * (ε_ref.cell[i,j] + ε_ref.cell[i-1,j]))
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

    @testset "m = 1 only (v0.6 guard)" begin
        @test_throws ErrorException solve_sim(gsolve, εsolve, λ; Nmodes=1, m=0,  report=false, plots=false)
        @test_throws ErrorException solve_sim(gsolve, εsolve, λ; Nmodes=1, m=-1, report=false, plots=false)
        @test_throws ErrorException adj_opt(design, [λ_target(λ)]; substrates, nClad, nCBG=nSiN, tCBG, m=0)
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
        @test se.driven.λ == BullseyeFDFD._select_mode(se.modes, λ; grid=gsolve,
                                                       inv_ϵr=se.inv_ϵr).λ   # dipole(grid) ⇒ selected resonance
        sp = solve_sim(gsolve, εsolve, λ; Nmodes=1, source=dipole(gsolve),
                       purcell=true, report=false, plots=false)
        @test isfinite(sp.driven.Purcell) && sp.driven.Purcell > 0   # LDOS bulk ratio
    end

    @testset "solve_sim — driven axis Eϕ is regular (m=1)" begin
        d = solve_sim(gsolve, εsolve, λ; source=dipole(gsolve),
                      report=false, plots=false, purcell=false).driven
        @test maximum(abs.(d.Eϕ[1, :] .+ 1im .* d.Eρ[1, :])) ≤ 1e-8 * maximum(abs, d.Eρ[1, :])
    end

    # Cavity modes keep ≲0.1 of their energy in the PML, PML/box modes ≳0.8. Among 8
    # modes on this coarse grid the max-CF one is a PML mode, which the old default drove.
    @testset "mode selection — PML-energy filter" begin
        s8 = solve_sim(gsolve, εsolve, λ; Nmodes=8, source=dipole(gsolve),
                       report=false, plots=false, purcell=false)
        pf = [md.pml_frac for md in s8.modes]
        @test all(0 .≤ pf .≤ 1)
        @test pf ≈ [BullseyeFDFD._pml_fraction(md, gsolve, s8.inv_ϵr) for md in s8.modes]
        @test pf[argmax(i -> s8.modes[i].Q, eachindex(s8.modes))] < 0.1   # fundamental
        @test maximum(pf) > 0.8                                          # a PML mode is present
        drv = only(filter(md -> md.λ == s8.driven.λ, s8.modes))
        @test drv.pml_frac ≤ 0.5 && drv.Q ≥ 30                           # never drives a PML mode
        @test BullseyeFDFD._select_index(s8.modes, λ, md -> md.CF, gsolve, s8.inv_ϵr) ==
              argmax(i -> s8.modes[i].CF, eachindex(s8.modes))           # explicit select overrides
        pml_only = filter(md -> md.pml_frac > 0.5, s8.modes)
        @test_logs (:warn, r"PML") BullseyeFDFD._select_index(pml_only, λ, nothing, gsolve, s8.inv_ϵr)
    end

    @testset "solve_sim — NA gates collection FOMs" begin
        s0 = solve_sim(gsolve, εsolve, λ;      Nmodes=1, report=false, plots=false)
        @test !hasproperty(s0.modes[1], :η)                  # NA omitted ⇒ no far-field crunch
        s1 = solve_sim(gsolve, εsolve, λ, 0.4; Nmodes=1, report=false, plots=false)
        @test hasproperty(s1.modes[1], :η) && isfinite(s1.modes[1].η)
    end

    # Near→far: Js = ẑ×H, Ms = −ẑ×E, Eθ = −(Lφ+Z₀Nθ), Eφ = Lθ−Z₀Nφ (e^{+jωt}, fields ∝ e^{−imφ}).
    # A brute-force Love integral over the plane calc_farfield samples (φ-sum without 2π/nφ)
    # is exactly (nφ/2)² × its intensity; the m=±1 mirror pair (Eφ → −Eφ) radiates identically.
    @testset "calc_farfield — Love transform & m=±1 mirror" begin
        sf = solve_sim(gsolve, εsolve, λ; Nmodes=1, report=false, plots=false)
        md, Ce = sf.modes[1], sf.Ce
        θs, I = calc_farfield(gsolve, md, Ce; m=1, resolution=7)

        _, Hρ, Hϕ, _ = calc_H_fields(gsolve, md, Ce)
        iz  = argmin(abs.(gsolve.zp .- BullseyeFDFD.default_monitors(gsolve).farfield))
        avg(v) = 0.5 .* (v[1:end-1] .+ v[2:end])
        nh  = min(gsolve.Nρ - 1, gsolve.idx_PMLρ)
        ρ   = gsolve.ρd[1:nh]; w = 0.5 .* (gsolve.ρp[2:nh+1].^2 .- gsolve.ρp[1:nh].^2)
        Eρs = md.Eρ[1:nh, iz];                              Eϕs = avg(md.Eϕ[:, iz])[1:nh]
        Hρs = avg(0.5 .* (Hρ[:, iz-1] .+ Hρ[:, iz]))[1:nh]; Hϕs = (0.5 .* (Hϕ[:, iz-1] .+ Hϕ[:, iz]))[1:nh]
        k0, Z0, nφ = 2π / md.λ, 376.73, 128
        love = map(θs) do θ
            Nx = Ny = Lx = Ly = 0.0im
            for φ in (0:nφ-1) .* (2π/nφ), i in 1:nh
                e = cis(-φ + k0 * ρ[i] * sin(θ) * cos(φ)) * w[i]
                Ex = Eρs[i]*cos(φ) - Eϕs[i]*sin(φ); Ey = Eρs[i]*sin(φ) + Eϕs[i]*cos(φ)
                Hx = Hρs[i]*cos(φ) - Hϕs[i]*sin(φ); Hy = Hρs[i]*sin(φ) + Hϕs[i]*cos(φ)
                Nx -= Hy*e; Ny += Hx*e; Lx += Ey*e; Ly -= Ex*e
            end
            abs2(Ly + Z0*cos(θ)*Nx) + abs2(cos(θ)*Lx - Z0*Ny)
        end
        @test love ≈ (nφ/2)^2 .* I  rtol=1e-8

        n   = gsolve.N
        mir = (; md..., m=-1, Eϕ=-md.Eϕ, raw_E=[md.raw_E[1:n]; -md.raw_E[n+1:2n]; md.raw_E[2n+1:3n]])
        @test calc_farfield(gsolve, mir, CeChOperators(gsolve, λ; m=-1)[1]; m=-1, resolution=7)[2] ≈ I  rtol=1e-10
    end

    # Uniform Eρ = Ez = 1 with Hϕ = h (H built through a stand-in Ce) gives Sz = Re(h)/2 and
    # Sρ = −Re(h)/2 everywhere, so the planes must return Sz·πR² and the wall Sρ·2πR·(z_top − z_bot)
    # exactly: the box faces tile the closed surface. The adjoint's copy must match calc_flux.
    @testset "calc_flux — box faces tile the surface; adjoint copy matches" begin
        n, h, Z0 = gsolve.N, 0.7 - 0.2im, 376.73
        raw_E = [ones(ComplexF64, n); zeros(ComplexF64, n); ones(ComplexF64, n)]
        f  = (; λ=1.0, raw_E, Eρ=reshape(raw_E[1:n], gsolve.Nρ, gsolve.Nz),
                Eϕ=reshape(raw_E[n+1:2n], gsolve.Nρ, gsolve.Nz), Ez=reshape(raw_E[2n+1:3n], gsolve.Nρ, gsolve.Nz))
        Ce = sparse(n .+ (1:n), 1:n, fill(-1im*Z0*2π*h, n), 3n, 3n)   # (Ce·E)ϕ = −iZ₀k₀h ⇒ Hϕ = h
        P_total, P_top, P_ρ, P_bot = calc_flux(gsolve, f, Ce; report=false)
        mon = BullseyeFDFD.default_monitors(gsolve)
        R   = gsolve.ρp[argmin(abs.(gsolve.ρd .- mon.right))]
        Hb  = gsolve.zp[argmin(abs.(gsolve.zp .- mon.top))] - gsolve.zp[argmin(abs.(gsolve.zp .- mon.bot))]
        @test P_top ≈ real(h)/2 * π * R^2  rtol=1e-12
        @test P_bot ≈ P_top                rtol=1e-12
        @test P_ρ   ≈ -real(h)/2 * 2π * R * Hb  rtol=1e-12

        sc = solve_sim(gsolve, εsolve, λ, 0.4; Nmodes=1, report=false, plots=false)
        md = sc.modes[1]
        @test BullseyeFDFD._flux_terms(md, gsolve, sc.Ce).top_frac ≈ md.P_top / md.P_total  rtol=1e-12
    end

    # -- adjoint ------------------------------------------------------------
    # _reciprocity_weight mirrors CeChOperators: W·Ch·Ce is complex-symmetric, so the left eigenvector
    # of A = ε⁻¹·Ch·Ce is W·ε·x. The residual catches drift if either side changes.
    @testset "CeChOperators ↔ _reciprocity_weight (left eigenvector)" begin
        for m in (1, -1)
            Ce, Ch = CeChOperators(gsolve, λ; m)
            iϵ, iμ = BullseyeFDFD.material_tensors(εsolve, gsolve)
            A  = iϵ * Ch * iμ * Ce
            md = first(BullseyeFDFD._solve_eig(A, gsolve, λ, 1, 1e-10, 3000))
            W  = BullseyeFDFD._reciprocity_weight(gsolve, λ; m)
            y  = BullseyeFDFD._left_eigvec(W, A, iϵ, md.raw_E, md.k_sq, gsolve)
            @test norm(A' * y .- conj(md.k_sq) .* y) / abs(md.k_sq) < 1e-8
        end
    end

    # Adjoint vs central finite differences on the frozen grid, for a width (eye) and a thickness (tCBG):
    # dk² (eigenvalue path), V (tracking box) and η (far field + flux, field path).
    @testset "adjoint gradients vs finite differences" begin
        B  = BullseyeFDFD; NA = 0.4
        gk = (; λPts=10, geoPts=6, PMLρ=0.3, PMLz_bot=0.3, PMLz_top=0.3, padρ=0.3, padz_bot=0.3, padz_top=0.3)
        p0 = [design; tCBG]; it = length(p0)
        g, M, εof, Ce, W = B._adjoint_setup(p0, substrates, nClad, nSiN, tCBG, false, λ, 1, gk;
                                            tspec=B._thickness_spec((; tCBG=(0.15, 0.22)), substrates))
        top(p) = g.cbg_bot + p[it]
        solve(p) = (ϵ = εof(p); iϵ = B.material_tensors(ϵ, g)[1]; A = iϵ * M;
                    (B._solve_eig(A, g, λ, 10, 1e-12, 3000), ϵ, iϵ, A))
        raw, ϵ, iϵ, A = solve(p0)
        md = B._select_mode(raw, λ; grid=g, inv_ϵr=iϵ)
        F  = B._shifted_lu(A, md.k_sq)
        dk2, y, b = B._dk2_dp(p0, md, g, A, iϵ, M, εof, W)
        ∇V, _ = B._dV_dp(p0, md, g, ϵ, iϵ, εof, A, F, y, b, dk2; top=top(p0), it)
        ∇η, _ = B._dη_dp(p0, md, g, iϵ, εof, A, F, y, b, dk2, Ce; NA)
        function obs(p)                                     # track the base mode by eigenvalue
            raw, ϵ, iϵ, _ = solve(p); m = argmin(r -> abs(r.k_sq - md.k_sq), raw)
            (m.k_sq, B._V_terms(m, g, ϵ, iϵ; top=top(p)).Veff, B.analyze_collection(m, g, Ce; NA).η)
        end
        for k in (1, it)
            δ = 1e-6; pp = copy(p0); pp[k] += δ; pm = copy(p0); pm[k] -= δ
            (kp, Vp, ηp), (km, Vm, ηm) = obs(pp), obs(pm)
            @test isapprox(dk2[k], (kp - km)/2δ; rtol=1e-4)
            @test isapprox(∇V[k],  (Vp - Vm)/2δ; rtol=1e-4)
            @test isapprox(∇η[k],  (ηp - ηm)/2δ; rtol=1e-4)
        end
    end
end
