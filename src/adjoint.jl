struct FOMTerm
    observable::Symbol
    kind::Symbol
    g::Function
    label::String
end

_lorentz_gate(o0, σ) = o -> (Δ = o - o0;      G = 1/(1 + (Δ/σ)^2); (G, -2Δ/σ^2 * G^2))
_sq_penalty(o0, w)   = o -> (Δ = (o - o0)/o0;                       (w*Δ^2, 2w*Δ/o0))

λ_target(λ0; σ=0.010) =
    FOMTerm(:λ, :gate, _lorentz_gate(λ0, σ), "λ_target($λ0; σ=$σ)")

Q_target(spec; w=1.0) =
    spec === :max ? FOMTerm(:Q, :reward,  o -> (o, 1.0),               "Q_target(:max)") :
                    FOMTerm(:Q, :penalty, _sq_penalty(float(spec), w), "Q_target($(float(spec)); w=$w)")

Fp_target(spec; w=1.0) =
    spec === :max ? FOMTerm(:Fp, :reward,  o -> (o, 1.0),               "Fp_target(:max)") :
                    FOMTerm(:Fp, :penalty, _sq_penalty(float(spec), w), "Fp_target($(float(spec)); w=$w)")

V_target(spec; w=1.0) =
    spec === :min ? FOMTerm(:V, :reward,  o -> (1/o, -1/o^2),           "V_target(:min)") :
                    FOMTerm(:V, :penalty, _sq_penalty(float(spec), w), "V_target($(float(spec)); w=$w)")

function evaluate_fom(terms, ovals)
    dJ_do = Dict{Symbol,Float64}()

    mult = [t for t in terms if t.kind !== :penalty]
    mv   = [t.g(ovals[t.observable]) for t in mult]
    for (i, t) in enumerate(mult)
        others = prod((first(mv[j]) for j in eachindex(mv) if j != i); init=1.0)
        dJ_do[t.observable] = get(dJ_do, t.observable, 0.0) - mv[i][2] * others
    end
    J = isempty(mv) ? 0.0 : -prod(first.(mv))

    for t in terms
        t.kind === :penalty || continue
        cost, dcost = t.g(ovals[t.observable])
        J += cost
        dJ_do[t.observable] = get(dJ_do, t.observable, 0.0) + dcost
    end
    return J, dJ_do
end

function _adjoint_setup(seed, substrates, nClad, nCBG, tCBG, slotted, λ, m, grid_kwargs)
    geom0  = build_geometry(seed, substrates, nClad, nCBG, tCBG; slotted)
    grid   = conformal_grid(geom0, λ; grid_kwargs...)
    Ce, Ch = CeChOperators(grid, λ; m=m)
    _, inv_μr = material_tensors(build_epsilon(geom0, grid), grid)
    M   = Ch * inv_μr * Ce
    εof = d -> build_epsilon(build_geometry(d, substrates, nClad, nCBG, tCBG; slotted), grid)
    return grid, M, εof
end

function _left_eigvec(A, k_sq, y0; steps=2)
    F = lu(A - (k_sq + 1e-6)*I)
    y = copy(y0)
    for _ in 1:steps
        y = F' \ y
        y ./= norm(y)
    end
    return y
end

function _zavgᵀ(d)
    Nρ, Nz = size(d); out = zeros(eltype(d), Nρ, Nz)
    @views out[:, 1] .+= d[:, 1]
    for j in 2:Nz
        @views out[:, j]   .+= 0.5 .* d[:, j]
        @views out[:, j-1] .+= 0.5 .* d[:, j]
    end
    return out
end

function _ρavgᵀ(d)
    Nρ, Nz = size(d); out = zeros(eltype(d), Nρ, Nz)
    @views out[1, :] .+= d[1, :]
    for i in 2:Nρ
        @views out[i, :]   .+= 0.5 .* d[i, :]
        @views out[i-1, :] .+= 0.5 .* d[i, :]
    end
    return out
end

function _dk2_deps_grid(dk2_dinv, inv_ϵr, grid)
    N = grid.N; d = diag(inv_ϵr)
    dρρ = reshape(dk2_dinv[1:N]     .* .-(d[1:N].^2),     grid.Nρ, grid.Nz)
    dϕϕ = reshape(dk2_dinv[N+1:2N]  .* .-(d[N+1:2N].^2),  grid.Nρ, grid.Nz)
    dzz = reshape(dk2_dinv[2N+1:3N] .* .-(d[2N+1:3N].^2), grid.Nρ, grid.Nz)
    return _zavgᵀ(dρρ) .+ _ρavgᵀ(dzz .+ _zavgᵀ(dϕϕ))
end

function _dk2_dparams(design, grad_ε, εof; δ=1e-7)
    grad = zeros(ComplexF64, length(design))
    for k in eachindex(design)
        dp = copy(design); dp[k] += δ
        dm = copy(design); dm[k] -= δ
        grad[k] = sum(grad_ε .* ((εof(dp) .- εof(dm)) ./ (2δ)))
    end
    return grad
end

function _dk2_dp(design, mode, grid, A, inv_ϵr, M, εof, y0; δ=1e-7)
    x = mode.raw_E
    y = _left_eigvec(A, mode.k_sq, y0)
    b = M * x
    dk2_dinv = conj(y) .* b ./ dot(y, x)
    grad_ε   = _dk2_deps_grid(dk2_dinv, inv_ϵr, grid)
    return _dk2_dparams(design, grad_ε, εof; δ=δ), y
end

function _fom_and_grad(k_sq, dk2, terms)
    k = sqrt(k_sq); kr = real(k); ki = imag(k); sgn = sign(ki)
    λ = 2π/kr; Q = kr/(2abs(ki))

    J, dJ_do = evaluate_fom(terms, Dict(:λ => λ, :Q => Q))   # add :V/:Fp values here later

    ∇ = zeros(Float64, length(dk2))
    for i in eachindex(dk2)
        dk = dk2[i]/(2k); dkr = real(dk); dki = imag(dk)
        dλ = -2π/kr^2 * dkr
        dQ = dkr/(2abs(ki)) - Q/abs(ki)*sgn*dki
        ∇[i] = get(dJ_do, :λ, 0.0)*dλ + get(dJ_do, :Q, 0.0)*dQ
    end
    # field observables later:  haskey(dJ_do,:V) && (∇ .+= dJ_do[:V] .* ∇V)
    return J, ∇, (; λ, Q)
end

function _select_mode(raw, λ_target, σ; Qfloor=30.0)
    cand = filter(m -> m.Q ≥ Qfloor, raw)
    isempty(cand) && return nothing
    score(m) = m.Q / (1 + ((m.λ - λ_target)/σ)^2)
    argmax(score, cand)
end

_fmt_obs(s, v) = s === :λ ? @sprintf("%.6f", v) :
                 s === :Q ? @sprintf("%.2f", v) :
                            @sprintf("%.6g", v)

function _open_log(save_data, seed, log_obs, meta)
    save_data === false && return nothing, nothing
    if save_data isa AbstractString && endswith(save_data, ".csv")
        path = save_data; mkpath(dirname(path))
    else
        base = save_data isa AbstractString ? save_data :
               (isempty(PROGRAM_FILE) ? pwd() : dirname(abspath(PROGRAM_FILE)))
        path = joinpath(mkpath(joinpath(base, "results")),
                        Dates.format(now(), "yymmdd-HHMMSS") * ".csv")
    end
    io = open(path, "w")
    gh = try readchomp(`git rev-parse --short HEAD`) catch; "unknown" end
    println(io, "# git=$gh")
    for line in meta; println(io, "# ", line); end
    println(io, "# seed=", seed)
    println(io, "eval,elapsed_s,J,", join(string.(log_obs), ","),
                ",", join(["p$i" for i in eachindex(seed)], ","))
    flush(io)
    return io, path
end

function _log_row(io, log_obs, n_eval, elapsed, J, ovals, x)
    io === nothing && return
    obs  = join((haskey(ovals, s) ? _fmt_obs(s, ovals[s]) : "NaN" for s in log_obs), ",")
    pstr = join((@sprintf("%.6f", xi) for xi in x), ",")
    @printf(io, "%d,%.1f,%.6f,%s,%s\n", n_eval, elapsed, J, obs, pstr)
    flush(io)
end

function adj_opt(seed_design, FOMs;
                 substrates, nClad, nCBG, tCBG, slotted=false,
                 λ=0.780, lower=nothing, upper=nothing,
                 Nmodes=5, m=1, grid_kwargs=(;),
                 σ_sel=0.020, Qfloor=30.0, save_data=false,
                 iterations=100, linesearch=BackTracking(),
                 options=Optim.Options(; iterations, outer_iterations=5,
                                         x_abstol=1e-3, f_reltol=1e-4,
                                         show_trace=true))

    grid, M, εof = _adjoint_setup(seed_design, substrates, nClad, nCBG, tCBG,
                                  slotted, λ, m, grid_kwargs)

    log_obs = unique([t.observable for t in FOMs])
    meta = ["cbg: nClad=$nClad  nCBG=$nCBG  tCBG=$tCBG",
            "substrates: " * join(["(n=$(s.n), h=$(s.height))" for s in substrates], "  "),
            "solve: grid_kwargs=$grid_kwargs  Nmodes=$Nmodes  m=$m",
            "FOM: " * join([t.label for t in FOMs], " · ")]
    logio, logpath = _open_log(save_data, seed_design, log_obs, meta)

    y_prev  = Ref{Union{Nothing,Vector{ComplexF64}}}(nothing)
    cache_x = Ref(Float64[]); cache = Ref{Any}(nothing); n_eval = Ref(0)
    t_start = time()

    function evaluate(x)
        if cache_x[] != x
            cache_x[] = copy(x); n_eval[] += 1
            ϵ      = εof(x)
            inv_ϵr = material_tensors(ϵ, grid)[1]
            A      = inv_ϵr * M
            raw = _solve_eig(A, grid, λ, Nmodes, 1e-6, 3000)
            rm  = _select_mode(raw, λ, σ_sel; Qfloor)
            if rm === nothing
                @printf("[%3d] no valid mode (no Q≥%.0f)  ‖Δp‖=%.3f\n",
                        n_eval[], Qfloor, norm(x .- seed_design))
                for (i, md) in enumerate(raw)
                    @printf("        raw %d: λ=%.1f nm  Q=%.1f\n", i, md.λ*1000, md.Q)
                end
                cache[] = (J = 1e3, ∇ = zeros(length(x)))
                _log_row(logio, log_obs, n_eval[], time()-t_start, 1e3, Dict{Symbol,Float64}(), x)
            else
                y0 = (y_prev[] !== nothing && length(y_prev[]) == length(rm.raw_E)) ?
                     y_prev[] : rm.raw_E
                dk2, y = _dk2_dp(x, rm, grid, A, inv_ϵr, M, εof, y0; δ=1e-7)
                J, ∇, ob = _fom_and_grad(rm.k_sq, dk2, FOMs)
                y_prev[] = y
                cache[]  = (; J, ∇)
                @printf("[%3d] J=%.4f  λ=%.2f nm  Q=%.1f\n", n_eval[], J, ob.λ*1000, ob.Q)
                _log_row(logio, log_obs, n_eval[], time()-t_start, J, Dict(:λ => ob.λ, :Q => ob.Q), x)
            end
        end
        return cache[]
    end
    f(x)     = evaluate(x).J
    g!(G, x) = (G .= evaluate(x).∇)

    result = try
        lower === nothing ?
            optimize(f, g!, seed_design, LBFGS(; linesearch), options) :
            optimize(f, g!, lower, upper, seed_design, Fminbox(LBFGS(; linesearch)), options)
    finally
        logio !== nothing && close(logio)
    end

    return (; design = Optim.minimizer(result), fom = Optim.minimum(result),
              grid, result, logfile = logpath)
end

