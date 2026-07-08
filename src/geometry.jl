function build_geometry(design, substrates, nClad, nCBG, tCBG; slotted=false)
    idx_eye = slotted ? 2 : 1
    N  = (length(design) - idx_eye - 1) ÷ 2

    @assert length(design) == idx_eye + 2N + 1 "length $(length(design)) doesn't fit a $(slotted ? "slotted" : "plain") bullseye" 

    cbg = (; eye    = (n = nCBG,  height=tCBG,  width = design[idx_eye]),
             trench = (n = nClad, height=0.0,   width = design[idx_eye+1 : 2 : idx_eye+2N-1]),
             ring   = (n = nCBG,  height=tCBG,  width = design[idx_eye+2 : 2 : idx_eye+2N]),
             buffer = (n = nClad, height=0.0,   width = design[end]))
    slotted && (cbg = (; slot = (n = nClad, height=0.0, width = design[1]), cbg...))

    k = length(substrates) - 1
    subnames = ntuple(i -> Symbol(:substrate, i), k)
    subs  = NamedTuple{subnames}(ntuple(i -> (n = substrates[i].n,
                                              height = substrates[i].height), k))
    wafer = (n = substrates[end].n, height = substrates[end].height)

    return (; cbg..., tCBG, stack = (; subs..., wafer))
end