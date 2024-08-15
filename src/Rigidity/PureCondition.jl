"""
    tiedown(g::Graphs.AbstractSimpleGraph, d::Integer; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}})

Preprocessing step to calculate pure condition of g. Transform g into a directed graph and add standard tiedown bars directed from vertices of the graph
to tiedown verts. If tiedown_verts are given, 1 bar is added to the first vert, 2 to the second, ... d to the dth.
"""
function tiedown(g::Graphs.AbstractSimpleGraph, d::Integer=2; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    @assert Graphs.nv(g) >= d "Number of vertices needs to be larger than dimension, but got $(Graphs.nv(g)) vertices and dimension $d."
    if isnothing(tiedown_verts)
        tiedown_verts = collect(1:d)
    else
        @assert tiedown_verts .<= Graphs.nv(g) "tiedown_verts need to consist of verts of g (1..$(nv(g)))"
        @assert tiedown_verts .>= 1 "tiedown_verts need to consist of verts of g (1..$(Graphs.nv(g)))"
        @assert length(tiedown_verts) == d "tiedown_verts needs to be vector of length d ($d), but got $(length(tiedown_verts))."
    end

    g_di = Graphs.SimpleDiGraph(g)

    # add tiedown bars as directed edges to graph
    for i in 1:d
        for j in 1:i
            Graphs.add_vertex!(g_di)
            Graphs.add_edge!(g_di, tiedown_verts[i], Graphs.nv(g_di))
        end
    end

    return g_di
end