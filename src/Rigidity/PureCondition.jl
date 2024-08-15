"""
    tiedown(g::Graphs.AbstractSimpleGraph, d::Integer; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}})

Preprocessing step to calculate pure condition of g. Transform g into a directed graph and add standard tiedown bars directed from vertices of the graph
to tiedown verts. If tiedown_verts are given, 1 bar is added to the dth vert, 2 to the (d-1)st, ... d to the 1st. If tiedown_verts is nothing, the first d vertices are chosen.

A directed graph is returned. This is because every edge from a vertex v to a vertex w corresponds to a non zero 1 × d submatrix of the rigidity matrix (exactly the submatrix indexed by the edge {v,w} and the vertex v).
This is, why we only add an edge  from a vertex v ̲~towards~ a tiedown vertex w as they add a 1 × d submatrix indexed with the edge {v,w} and the vertex v."""
function tiedown(g::Graphs.AbstractSimpleGraph, d::Integer=2; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    @assert Graphs.nv(g) >= d "Number of vertices needs to be larger than dimension, but got $(Graphs.nv(g)) vertices and dimension $d."
    if isnothing(tiedown_verts)
        tiedown_verts = collect(1:d)
    else
        @assert all(tiedown_verts .<= Graphs.nv(g)) "tiedown_verts need to consist of verts of g (1..$(nv(g)))"
        @assert all(tiedown_verts .>= 1) "tiedown_verts need to consist of verts of g (1..$(Graphs.nv(g)))"
        @assert length(tiedown_verts) == d "tiedown_verts needs to be vector of length d ($d), but got $(length(tiedown_verts))."
    end

    g_di = Graphs.SimpleDiGraph(g)

    # add tiedown bars as directed edges to graph
    for i in 1:d
        for j in 1:(d+1)-i
            Graphs.add_vertex!(g_di)
            Graphs.add_edge!(g_di, tiedown_verts[i], Graphs.nv(g_di))
        end
    end

    return g_di
end

function tiedown(poly::AbstractEmbOrCombPolyhedron, d::Integer=3; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    if !isnothing(tiedown_verts)
        return tiedown(Graphs.SimpleGraph(poly), d; tiedown_verts=tiedown_verts)
    end

    if d < max(length(get_facets(poly)))
        i = findfirst(f -> length(f) >= d, get_facets(poly))
        f = get_facets(poly)[i]
        tiedown_verts = f[1:d]
        return tiedown(Graphs.SimpleGraph(poly), d; tiedown_verts=tiedown_verts)
    end

    return tiedown(Graphs.SimpleGraph(poly), d)
end

# tiedown factor of the standard tiedown at the vertices tiedown_verts as in White and Whiteley 1983. Tiedown is applied like in the function tiedown!
#  Is only true, if induced subgraph by tiedown_verts is d-isostatic
function tiedown_factor(g::Graphs.AbstractSimpleGraph, d::Integer=2; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    @assert is_isostatic(g[tiedown_verts], d)

    @assert Graphs.nv(g) >= d "Number of vertices needs to be larger than dimension, but got $(Graphs.nv(g)) vertices and dimension $d."
    if isnothing(tiedown_verts)
        tiedown_verts = collect(1:d)
    else
        @assert all(tiedown_verts .<= Graphs.nv(g)) "tiedown_verts need to consist of verts of g (1..$(nv(g)))"
        @assert all(tiedown_verts .>= 1) "tiedown_verts need to consist of verts of g (1..$(Graphs.nv(g)))"
        @assert length(tiedown_verts) == d "tiedown_verts needs to be vector of length d ($d), but got $(length(tiedown_verts))."
    end

    B = BracketAlgebra(d, Graphs.nv(g) + binomial(d + 1, 2))

    factors = Vector{Int}[]
    next = Graphs.nv(g) + 1

    for i in eachindex(tiedown_verts)
        factor = tiedown_verts[1:i]
        for _ in i+1:d+1
            push!(factor, next)
            next += 1
        end
        push!(factors, factor)
    end

    tabloids = map(factor -> Tabloid([factor]), factors)
    return prod(bracket_monomial(tabloid, B) for tabloid in tabloids)
end

# reduce the graph g by removing all edges in edges and their reverse edges and removing all outward edges from v
# this corresponds to deleting the rows of the rigidity matrix indexed by edges and the columns indexed by v as during Laplace expansion.
function reduction!(g::Graphs.SimpleDiGraph, v::Integer, edges::Vector{<:Graphs.AbstractEdge}, d::Integer)
    @assert all(map(e -> Graphs.has_edge(g, Graphs.src(e), Graphs.dst(e)), edges)) "edges needs to be subset of edges of g."
    @assert length(edges) == d "Number of edges needs to be equal to dimension, but got $(length(edges)) edges and dimension $d."
    @assert all(Graphs.src.(edges) .== v) "All edges need to start at vertex $v, but edge sources are $(src.(edges))."

    remove = union(edges, reverse.(edges), [Graphs.Edge(v, w) for w in Graphs.outneighbors(g, v)])
    for e in remove
        Graphs.rem_edge!(g, e)
    end

    return g
end

function reduction(g::Graphs.SimpleDiGraph, v::Integer, edges::Vector{<:Graphs.AbstractEdge}, d::Integer)
    return reduction!(deepcopy(g), v, edges, d)
end

# g represents a submatrix of the rigidity matrix of some original graph with vertices 1:nv(g) : 
# - If a vertex v in g still has outedges, this means that the columns corresponding to v have not been eliminated via Laplace expansion in a previous step.
# - If an edge (v,w) in g still exists, this means that the row corresponding to {v,w} has not been eliminated via Laplace expansion in a previous step.
# This function computes the sign of the bracket expression that occurs when deleting the rows indexed by edges and the columns indexed by v
# from the submatrix of the original rigidity matrix represented by the graph g.
function sign(g::Graphs.SimpleDiGraph, v::Integer, edges::Vector{<:Graphs.AbstractEdge}, d::Integer)
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

# recursively calculate the condition for the directed graph g to be infinitesimally flexible as an expression in the bracket algebra B.
# This is done by iterating Laplace expansion of the rigidity matrix of g by expanding aling the columns corresponding to the vertices.
# The expansion can be read off the graph without actually constructing the rigidity matrix.
# For reference see https://omni.wikiwand.com/en/articles/Laplace_expansion#General_statement 
function condition(g::Graphs.SimpleDiGraph, B::BracketAlgebra)
    d = B.d

    if Graphs.ne(g) == 0
        return 1
    end

    relevant_verts = filter(v -> length(Graphs.outneighbors(g, v)) > 0, Graphs.vertices(g))

    # select the vertex with smallest defect: d - (outdegree - indegree).
    (defect, i) = findmin(v -> d - (length(Graphs.outneighbors(g, v)) - length(Graphs.inneighbors(g, v))), relevant_verts)
    v = relevant_verts[i]

    if defect <= 0
        # if defect == 0 that means the vertex with minimum defect satisfies (outdegree - indegree) = d. 
        # This means the rigidity matrix can be rearranged as a upper left triangular block matrix with a d×d block in the upper left corner. 
        # Thus, we get the determinant of the matrix by multiplying the determinant of this block with the determinant of the matrix after deleting rows and columns of the block.
        # The determinant of the block is a bracket expression with v and d outwards neighbors of v that are not inwards neighbors of v.

        # edges from v that don't have a reverse edge 
        edges = setdiff([Graphs.Edge(v, w) for w in Graphs.outneighbors(g, v)], [Graphs.Edge(v, w) for w in Graphs.inneighbors(g, v)])[1:d]

        # recursive call. The determinant of the rigiditymatrix is via Laplace: ± [v, e1_2, …, ed_2]  * (determinant of matrix after deleting rows corresponding to edges and columns corresponding to v)
        return sign(g, v, edges, d) * bracket_monomial(Tabloid([pushfirst!([Graphs.dst(e) for e in edges], v)]), B) * condition(reduction(g, v, edges, d), B)
    elseif defect > 0
        # If defect > 0 the determinant has to be calculated using Laplace expansion that involve more than one nonzero summand. 
        # See Wikipedia article.

        edges_onlyout = setdiff([Graphs.Edge(v, w) for w in Graphs.outneighbors(g, v)], [Graphs.Edge(v, w) for w in Graphs.inneighbors(g, v)])
        edges_both = [Graphs.Edge(v, w) for w in intersect(Graphs.outneighbors(g, v), Graphs.inneighbors(g, v))]

        # recursive call for Laplace expansion
        sum_index = map(edges -> union(edges, edges_onlyout), combinations(edges_both, defect))
        return sum(edges -> sign(g, v, edges, d) * bracket_monomial(Tabloid([pushfirst!([Graphs.dst(e) for e in edges], v)]), B) * condition(reduction(g, v, edges, d), B), sum_index)
    elseif defect < 0
        # If defect < 0 the determinant is zero as after Laplace expansion along the columns of v and d outedges of v, the matrix has a zero row.
        return 0
    end
end

function condition(g::Graphs.SimpleDiGraph, d::Integer=2)
    B = BracketAlgebra(g, d)

    return condition(g, B)
    # return Groebner.normalform(reduced_groebner_basis!(B), condition(g, B), ordering=B.ordering)
end

function condition(g::Graphs.AbstractSimpleGraph, d::Integer=2; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    if !isnothing(tiedown_verts) && is_isostatic(g[tiedown_verts], d)
        tiedown_fac = tiedown_factor(g, d; tiedown_verts=tiedown_verts)
        return condition(tiedown(g, d; tiedown_verts=tiedown_verts), d) / tiedown_fac
    end

    return condition(tiedown(g, d; tiedown_verts=tiedown_verts), d)
end

function condition(poly::AbstractEmbOrCombPolyhedron, d::Integer=3; tiedown_verts::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    g = Graphs.SimpleGraph(poly)

    if !isnothing(tiedown_verts)
        return condition(g, d; tiedown_verts=tiedown_verts)
    end

    if d < max(length(get_facets(poly)))
        i = findfirst(f -> length(f) >= d, get_facets(poly))
        f = get_facets(poly)[i]
        tiedown_verts = f[1:d]
    else
        tiedown_verts = collect(1:d)
    end

    return condition(g, d, tiedown_verts=tiedown_verts)
end