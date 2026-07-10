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
design = [0.355000,             # eye
          0.100000, 0.370000,   # trench 1, ring 1
          0.100000, 0.370000,   # trench 2, ring 2
          0.100000, 0.370000,   # trench 3, ring 3
          0.100000, 0.370000,   # trench 4, ring 4
          0.100000, 0.370000,   # trench 5, ring 5
          0.100000, 0.370000,   # trench 6, ring 6
          0.100000, 0.370000,   # trench 7, ring 7
          0.100000, 0.370000,   # trench 8, ring 8
          0.100000, 0.370000,   # trench 9, ring 9
          0.100000, 0.370000,   # trench 10, ring 10
          0.100000];            # buffer

geometry = build_geometry(design, substrates, nClad, nSiN, tCBG);

# --- mesh ---
λ_target = 0.780;
grid = conformal_grid(geometry, λ_target);

# --- bullseye generation ---
ϵ = build_epsilon(geometry, grid);
show_sim(ϵ, grid);

# --- simulation setup ---
@time sim = solve_sim(grid, ϵ, λ_target;
                      Nmodes=3);

NA = 0.4;
@time sim = solve_sim(grid, ϵ, λ_target, NA;
                      Nmodes=3);

λ_drive = 0.781;
src = dipole(grid, λ_drive);

sim = solve_sim(grid, ϵ, λ_target;
                source=src);

src = dipole(grid);

@time sim = solve_sim(grid, ϵ, λ_target, NA;
                      Nmodes=3, source=src);