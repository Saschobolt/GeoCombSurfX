# ################## Simple combinatorial Graphs
# abstract type AbstractSimpleGraph{T<:Integer} end

# mutable struct Graph{T<:Integer} <: AbstractSimpleGraph{T}
#     verts::Vector{T}
#     edges::Vector{Vector{T}}

#     """
#     Graph(verts::Vector{<:Integer}, edges::Vector{<:Vector{<:Integer}})

#     TBW
#     """
#     function Graph(verts::Vector{<:Integer}, edges::Vector{<:Vector{<:Integer}})
#         @assert all([length(e) == 2 for e in edges]) "Edges need to be of length 2."
#         @assert unique(Set.(edges)) == Set.(edges) "There are duplicate edges."
#         @assert issubset(Set(vcat(edges...)), Set(verts)) "Edges need to consist of vertices of the graph."
#         @assert Set(verts) == Set(1:length(verts)) "Verts have to be 1,..., n."

#         T = typeof(verts[1])
#         return new{T}(sort(verts), edges)
#     end

#     function Graph(;verts = nothing, edges = nothing)
#         if isnothing(edges)
#             edges = Vector{Int64}[]
#         end
#         if isnothing(verts)
#             if edges == []
#                 verts = Int64[]
#             else
#                 verts = collect(1:max(vcat(edges...)...))
#             end
#         end

#         return Graph(verts, edges)
#     end
# end


# """
#     get_verts(g::AbstractSimpleGraph)

# TBW
# """
# get_verts(g::AbstractSimpleGraph) = deepcopy(g.verts)

# function set_verts!(g::Graph, verts::Vector{<:Integer})
#     g.verts = verts
# end

# """
#     get_edges(g::AbstractSimpleGraph)

# TBW
# """
# get_edges(g::AbstractSimpleGraph) = deepcopy(g.edges)

# """
#     set_edges!(g::AbstractSimpleGraph, edges::Vector{<:Vector{<:Integer}})

# TBW
# """
# function set_edges!(g::AbstractSimpleGraph, edges::Vector{<:Vector{<:Integer}})
#     if setdiff(vcat(edges...), get_verts(g)) != []
#         error("Edges need to consist of vertices of the graph.")
#     end
#     g.edges = edges
# end

# """
#     no_concomponents(g::AbstractSimpleGraph)

# TBW
# """
# function no_concomponents(g::AbstractSimpleGraph)
#     verts = collect(1:max(vcat(get_edges(g)...)...))
#     sol = 0

#     visited = [false for v in verts]

#     function update!(k)
#         visited[k] = true

#         inc_edges = filter(e -> k in e, get_edges(g))
#         if inc_edges == []
#             return
#         end

#         adj = unique(filter(l -> !visited[l], union(inc_edges...)))
#         for l in adj
#             update!(l)
#         end
#     end

#     for k in verts
#         if !visited[k]
#             sol = sol + 1
#             update!(k)
#         end
#     end

#     return sol
# end


################## Embedded Graphs

mutable struct Framework{S<:Real,T<:Integer} <: AbstractEmbeddedGraph{S,T}
    verts::Matrix{S}
    edges::Vector{MVector{2,T}}

    """
    Framework(verts::Matrix{<:Real}, edges::Vector{<:Vector{<:Integer}})

    TBW
    """
    function Framework(verts::AbstractMatrix{<:Real}, edges::AbstractVector{<:AbstractVector{<:Integer}})
        # combgraph = Graph(edges=edges)
        # if no_concomponents(combgraph) != 1
        #     error("Graph is not connected.")
        # end
        if size(verts)[2] < max(vcat(Vector.(edges)...)...)
            error("Not every vertex of an edge has a coordinate assigned to it.")
        end

        S = typeof(promote(verts[1, 1], 1.0)[1])
        T = typeof(edges[1][1])
        return new{S,T}(verts, edges)
    end

    function Framework(; verts::Union{AbstractMatrix{<:Real},Nothing}, edges::Union{AbstractVector{AbstractVector{<:Integer}},Nothing})
        if isnothing(edges)
            edges = Vector{Int64}[]
        end
        if isnothing(verts)
            verts = rand(Float64, (3, max(vcat(Vector.(edges)...)...)))
        end

        return Framework(verts, edges)
    end
end

get_verts(f::AbstractEmbeddedGraph) = deepcopy(f.verts)


get_edges(f::AbstractEmbeddedGraph) = deepcopy(f.edges)

function adjacency_matrix(f::AbstractEmbeddedGraph)
    m = zeros(Int, size(get_verts(f))[2], size(get_verts(f))[2])
    for e in get_edges(f)
        m[e[1], e[2]] = 1
        m[e[2], e[1]] = 1
    end
    return m
end

Graphs.SimpleGraph(g::AbstractEmbeddedGraph) = Graphs.SimpleGraph(adjacency_matrix(g))

convert(::Type{T}, f::AbstractEmbeddedGraph) where {T<:Graphs.AbstractSimpleGraph} = T(adjacency_matrix(f))

get_edges(g::Graphs.AbstractSimpleGraph) = [[e.src, e.dst] for e in Graphs.edges(g)]

get_verts(g::Graphs.AbstractSimpleGraph) = collect(1:Graphs.nv(g))

function Framework(g::Graphs.AbstractSimpleGraph, d::Integer=3)
    verts = rand(Float64, (d, size(Graphs.adjacency_matrix(g))[1]))
    return Framework(verts, get_edges(g))
end

"""
    set_verts!(g::AbstractEmbeddedGraph, verts::Matrix{<:Real})

TBW
"""
function set_verts!(g::AbstractEmbeddedGraph, verts::AbstractMatrix{<:Real})
    @assert size(verts)[2] >= max(vcat(get_edges(g)...)...) "Not every vertex of an edge has a coordinate assigned to it."
    g.verts = verts
end

"""
    dimension(g::AbstractEmbeddedGraph)

TBW
"""
dimension(g::AbstractEmbeddedGraph) = size(get_verts(g))[1]

function display(g::Framework)
    print(
        """$(typeof(g)) in $(dimension(g))-space with $(size(get_verts(g))[2]) vertices and $(length(get_edges(g))) edges.
        Edges: $(get_edges(g))
        """
    )
end

function edge_lengths(f::Framework)
    [sqrt(sqdist(get_verts(f)[:, e[1]], get_verts(f)[:, e[2]])) for e in get_edges(f)]
end

function henneberg_extension!(g::Graphs.AbstractSimpleGraph, verts::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}}; check::Bool=true)
    if check
        @assert all([length(e) == 2 for e in edges]) "Edges need to be of length 2."
        @assert issubset(Set.(edges), Set.(get_edges(g))) "edges need to consist of edges of g."
        @assert issubset(verts, get_verts(g)) "verts need to consist of vertices of g."
    end

    @assert all([v in verts for v in vcat(edges...)]) "edges need to connect vertices in verts."

    for e in edges
        Graphs.rem_edge!(g, e[1], e[2])
    end

    Graphs.add_vertex!(g)

    for v in verts
        Graphs.add_edge!(g, v, Graphs.nv(g))
    end

    return g
end

"""
    henneberg_extension!(f::Framework, verts::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}}, p::Union{AbstractVector{<:Real},Nothing}=nothing; check::Bool=true)

Perform a Henneberg extension on the framework f by adding a new vertex with coordinates p, connecting it to the vertices verts and removing the edges in edges.
"""
function henneberg_extension!(f::Framework, verts::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}}, p::Union{AbstractVector{<:Real},Nothing}=nothing; check::Bool=true)
    g = SimpleGraph(f)
    henneberg_extension!(g, verts, edges, check=check)

    if isnothing(p)
        f.verts = hcat(f.verts, rand(Float64, size(f.verts)[1]))
    else
        f.verts = hcat(f.verts, p)
    end


    f.edges = get_edges(Framework(g))

    return f
end

function henneberg_extension(f::Framework, verts::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}}, p::Union{AbstractVector{<:Real},Nothing}=nothing; check::Bool=true)
    return henneberg_extension!(deepcopy(f), verts, edges, p, check=check)
end