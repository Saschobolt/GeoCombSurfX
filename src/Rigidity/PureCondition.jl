"""
    tiedown(g::Graphs.AbstractSimpleGraph, d::Integer; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}})

Preprocessing step to calculate pure condition of g. Transform g into a directed graph and add standard tiedown bars directed from vertices of the graph
to tiedown verts. If tiedown_verts are given, 1 bar is added to the first vert, 2 to the second, ... d to the dth.

A directed graph is returned. This is because every edge from a vertex v to a vertex w corresponds to a non zero 1 × d submatrix of the rigidity matrix (exactly the submatrix indexed by the edge {v,w} and the vertex v).
This is, why we only add an edge  from a vertex v ̲~towards~ a tiedown vertex w as they add a 1 × d submatrix indexed with the edge {v,w} and the vertex v."""
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

# reduce the graph g by removing all edges in edges and their reverse edges and removing all outward edges from v
# this corresponds to deleting the rows of the rigidity matrix indexed by edges and the columns indexed by v as during Laplace expansion.
function reduction!(g::Graphs.SimpleDiGraph, v::Integer, edges::Vector{<:Graphs.AbstractEdge}, d::Integer)
    @assert all(map(e -> Graphs.has_edge(g, Graphs.src(e), Graphs.dst(e)), edges)) "edges needs to be subset of edges of g."
    @assert length(edges) == d "Number of edges needs to be equal to dimension, but got $(length(e)) edges and dimension $d."
    @assert all(Graphs.src.(edges) .== v) "All edges need to start at vertex $v, but edge sources are $(src.(e))."

    remove = union(e, reverse.(e), [Graphs.Edge(v, w) for w in Graphs.outneighbors(g, v)])
    for e in remove
        Graphs.rem_edge!(g, e)
    end

    return g
end

# g represents a submatrix of the rigidity matrix of some original graph with vertices 1:nv(g) : 
# - If a vertex v in g still has outedges, this means that the columns corresponding to v have not been eliminated via Laplace expansion in a previous step.
# - If an edge (v,w) in g still exists, this means that the row corresponding to {v,w} has not been eliminated via Laplace expansion in a previous step.
# This function computes the sign of the bracket expression that occurs when deleting the rows indexed by edges and the columns indexed by v
# from the submatrix of the original rigidity matrix represented by the graph g.
function _sign(g::Graphs.SimpleDiGraph, v::Integer, edges::Vector{<:Graphs.AbstractEdge}, d::Integer)
    @assert all(map(e -> Graphs.has_edge(g, Graphs.src(e), Graphs.dst(e)), edges)) "edges needs to be subset of edges of g."
    @assert length(edges) == d "Number of edges needs to be equal to dimension, but got $(length(e)) edges and dimension $d."
    @assert all(Graphs.src.(edges) .== v) "All edges need to start at vertex $v, but edge sources are $(src.(e))."

    # vertices that still have columns in the rigidity matrix are those that have outgoing edges
    relevant_verts = filter(v -> length(Graphs.outneighbors(g, v)) > 0, Graphs.vertices(g))
    i = indexin(v, relevant_verts)[1]
    # indices of the columns indexed by v in the remaining rigidity matrix
    col_indices = d*(i-1)+1:d*(i-1)+d

    undirected_edges = unique(filter(e -> Graphs.src(e) < Graphs.dst(e), union(Graphs.edges(g), reverse.(Graphs.edges(g)))))

    # indices of the rows indexed by edges in the remaining rigidity matrix. Wlog we assume that the rows of the rigidity matrix are sorted lexicographically.
    sort!(undirected_edges, by=e -> (Graphs.src(e), Graphs.dst(e)))
    row_indices = indexin(map(e -> (Graphs.src(e) > Graphs.dst(e)) ? Graphs.reverse(e) : e, edges), undirected_edges)

    return (-1)^(sum(row_indices) + sum(col_indices))
end

