using Revise; using BullseyeFDFD;

# --- materials ---
nClad = 1.00;
nSiN  = 2.01066;
nSiO2 = 1.45375;
nSi   = 3.69476 - 1im*0.00482;

# --- layer heights (μm) ---
tCBG  = 0.17886;
tSiO2 = 0.75057;

# --- substrates ---
substrates = [(n = nSiO2,   height = tSiO2), 
              (n = nSi,     height = Inf)];

# --- design ---
design = [0.037500, 0.273280,   # slot, eye
          0.099090, 0.378089,   # trench 1, ring 1
          0.147451, 0.357824,   # trench 2, ring 2
          0.162671, 0.386911,   # trench 3, ring 3
          0.158265, 0.388490,   # trench 4, ring 4
          0.169946, 0.362279,   # trench 5, ring 5
          0.181444, 0.347768,   # trench 6, ring 6
          0.189631];            # buffer

geometry = build_geometry(design, substrates, nClad, nSiN, tCBG, slotted=true);

λ_target = 0.7825;
grid = conformal_grid(geometry, λ_target);

ϵ = build_epsilon(geometry, grid);

show_sim(ϵ, grid);

src = dipole(grid);

NA = 0.4;

@time sim = solve_sim(grid, ϵ, λ_target, NA;
                      Nmodes=3, source=src);