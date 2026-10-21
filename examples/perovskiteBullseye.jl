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
design = [0.037500, 0.300000,   # slot, eye
          0.100000, 0.380000,   # trench 1, ring 1
          0.100000, 0.380000,   # trench 2, ring 2
          0.100000, 0.380000,   # trench 3, ring 3
          0.100000, 0.380000,   # trench 4, ring 4
          0.100000, 0.380000,   # trench 5, ring 5
          0.100000, 0.380000,   # trench 6, ring 6
          0.100000];            # buffer

geometry = build_geometry(design, substrates, nClad, nSiN, tCBG, slotted=true);

λ0 = 0.7825;
grid = conformal_grid(geometry, λ0);

ϵ = build_epsilon(geometry, grid);

show_sim(ϵ, grid);

NA = 0.4;

@time sim = solve_sim(grid, ϵ, λ0, NA;
                      Nmodes=3);