# Rewrite the generated blocks of README.md from saved data. Standard library only; safe to re-run.
#
#   julia assets/update_readme.jl            # rewrite README.md
#   julia assets/update_readme.jl --check    # only report blocks that are out of date (exit 1 if any)
#
# A generated block is
#     <!-- readme:PATH -->
#     …
#     <!-- /readme:PATH -->
# where PATH is a data file relative to the repo root: a .txt is inserted as a fenced verbatim block,
# a .csv as a Markdown table (headers and number formats from HEADERS/FORMATS; a column with one value
# in every row moves to a line under the table). The data carries its provenance (a "# provenance:
# git=… src_clean=…" first line in a .txt; git/src_clean/… columns in a .csv). It is left out of the
# README, and the script warns when the data came from a commit whose src/ or Project.toml differs from HEAD.
using Printf

const ROOT      = dirname(@__DIR__)
const README    = joinpath(ROOT, "README.md")
const PROV_COLS = ("git", "src_clean", "version", "date")
const MARKER    = r"<!-- readme:(\S+) -->\n.*?<!-- /readme:\1 -->"s

# README header of a CSV column; unlisted columns keep their CSV name
const HEADERS = Dict(
    # validation/comsol/results/summary.csv
    "sweep" => "Sweep", "quantity" => "Quantity", "same_mode" => "Same-mode points",
    "ratio_mean" => "FDFD/COMSOL mean", "ratio_std" => "FDFD/COMSOL std",
    "max_abs_resid_pct" => "Max. difference (%)", "r" => "Correlation r",
    # validation/adjoint/results/summary.csv
    "observable" => "Observable", "params" => "Parameters", "δ_fd" => "FD step δ",
    "median_relerr" => "Median rel. error", "max_relerr" => "Max. rel. error")

# Printf format of a CSV column's numbers; unlisted columns print integers as-is, other numbers as %.4g
const FORMATS = Dict(
    "ratio_mean" => "%.4f", "ratio_std" => "%.4f", "max_abs_resid_pct" => "%.2f", "r" => "%.4f",
    "δ_fd" => "%.0e", "median_relerr" => "%.1e", "max_relerr" => "%.1e")

git_ok(args) = success(pipeline(`git -C $ROOT $args`; stdout=devnull, stderr=devnull))

function check_provenance(path, git, src_clean)
    if git === nothing
        @warn "$path: no provenance recorded"
    elseif !git_ok(`cat-file -e $(git * "^{commit}")`)
        @warn "$path: produced at unknown commit $git"
    elseif !git_ok(`diff --quiet $git HEAD -- src Project.toml`)
        changed = readlines(`git -C $ROOT diff --name-only $git HEAD -- src Project.toml`)
        @warn "$path: produced at $(git[1:7]); since then $(join(changed, ", ")) changed — re-run its producer if that affects it"
    end
    src_clean == "false" && @warn "$path: produced with uncommitted changes under src/ or Project.toml"
end

function render_txt(path)
    lines = readlines(joinpath(ROOT, path))
    m = isempty(lines) ? nothing : match(r"^# provenance: git=(\S+) src_clean=(\S+)", lines[1])
    check_provenance(path, m === nothing ? nothing : m[1], m === nothing ? nothing : m[2])
    m === nothing || popfirst!(lines)
    return "```\n" * join(lines, "\n") * "\n```"
end

header(name) = get(HEADERS, name, name)

function fmt(name, s)
    x = tryparse(Float64, s)
    x === nothing || occursin(r"^-?\d+$", s) ? s :
    haskey(FORMATS, name) ? Printf.format(Printf.Format(FORMATS[name]), x) : @sprintf("%.4g", x)
end

function render_csv(path)
    text = read(joinpath(ROOT, path), String)
    occursin('"', text) && error("$path: quoted CSV fields are not supported")
    rows = split.(filter(!isempty, split(text, '\n')), ',')
    head = rows[1]; body = rows[2:end]
    col(name) = (i = findfirst(==(name), head); i === nothing ? nothing : unique(r[i] for r in body))
    gits = col("git"); clean = col("src_clean")
    gits === nothing ? check_provenance(path, nothing, nothing) :
        foreach(g -> check_provenance(path, g, clean === nothing ? nothing : ("false" in clean ? "false" : "true")), gits)
    keep  = findall(h -> !(h in PROV_COLS), head)
    fixed = length(body) > 1 ? filter(i -> allequal(r[i] for r in body), keep) : Int[]
    cols  = setdiff(keep, fixed)
    tbl = ["| " * join(header.(head[cols]), " | ") * " |", "|" * repeat(" --- |", length(cols))]
    append!(tbl, ["| " * join(fmt.(head[cols], r[cols]), " | ") * " |" for r in body])
    isempty(fixed) ||
        push!(tbl, "", "*" * join(["$(header(head[i])): $(fmt(head[i], body[1][i]))" for i in fixed], " · ") * "*")
    return join(tbl, "\n")
end

function render(path)
    isfile(joinpath(ROOT, path)) || error("README block $path: no such file (run its producer first)")
    endswith(path, ".txt") ? render_txt(path) :
    endswith(path, ".csv") ? render_csv(path) : error("$path: only .txt and .csv blocks are supported")
end

old = read(README, String)
new = replace(old, MARKER => s -> begin
    path = match(MARKER, s)[1]
    "<!-- readme:$path -->\n" * render(path) * "\n<!-- /readme:$path -->"
end)

nblocks = length(collect(eachmatch(MARKER, old)))
if "--check" in ARGS
    new == old ? println("README.md: $nblocks generated blocks up to date") :
                 (println("README.md: generated blocks are out of date — run julia assets/update_readme.jl"); exit(1))
elseif new == old
    println("README.md: $nblocks generated blocks already up to date")
else
    write(README, new)
    println("README.md: rewrote $nblocks generated blocks")
end
