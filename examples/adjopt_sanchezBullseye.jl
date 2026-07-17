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
res  = adj_opt(design, FOMs; substrates, nClad, nCBG=nSiN, tCBG, λ=λ0,
               lower=fill(0.05, length(design)), upper=fill(1.00, length(design)),
               Nmodes=5, iterations=50, grid_kwargs = (; λPts=20, geoPts=8));
