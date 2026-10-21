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

λ0 = 0.780;

FOMs = [λ_target(λ0), Q_target(:max)];

grid_kwargs = (; λPts=24, geoPts=8);
opt_kwargs  = (; nClad, nCBG=nSiN, λ=λ0,
                 lower=fill(0.05, length(design)), upper=fill(1.00, length(design)),
                 thickness=(; tCBG=(0.15, 0.22), substrate1=(0.65, 0.85)),
                 Nmodes=10, iterations=50, grid_kwargs,
                 save_data=@__DIR__);

res1 = adj_opt(design,      FOMs; substrates,      tCBG,      opt_kwargs...);
res2 = adj_opt(res1.design, FOMs; res1.substrates, res1.tCBG, opt_kwargs...);
res3 = adj_opt(res2.design, FOMs; res2.substrates, res2.tCBG, opt_kwargs...);

geometry = build_geometry(res3.design, res3.substrates, nClad, nSiN, res3.tCBG);
grid     = conformal_grid(geometry, λ0; λPts=36, geoPts=24);
ϵ        = build_epsilon(geometry, grid);
sim      = solve_sim(grid, ϵ, λ0; Nmodes=10);
