import Base.Multimedia.display

################## Simple combinatorial Graphs
abstract type AbstractSimpleGraph{T<:Integer} end

mutable struct Graph{T<:Integer} <: AbstractSimpleGraph{T}
    verts::Vector{T}
    edges::Vector{Vector{T}}

    """
    Graph(verts::Vector{<:Integer}, edges::Vector{<:Vector{<:Integer}})

    TBW
    """
    function Graph(verts::Vector{<:Integer}, edges::Vector{<:Vector{<:Integer}})
        @assert all([length(e) == 2 for e in edges]) "Edges need to be of length 2."
        @assert unique(Set.(edges)) == Set.(edges) "There are duplicate edges."
        @assert issubset(Set(vcat(edges...)), Set(verts)) "Edges need to consist of vertices of the graph."
        @assert Set(verts) == Set(1:length(verts)) "Verts have to be 1,..., n."

        T = typeof(verts[1])
        return new{T}(sort(verts), edges)
    end
    
    function Graph(;verts = nothing, edges = nothing)
        if isnothing(edges)
            edges = Vector{Int64}[]
        end
        if isnothing(verts)
            if edges == []
                verts = Int64[]
            else
                verts = collect(1:max(vcat(edges...)...))
            end
        end
    
        return Graph(verts, edges)
    end
end


"""
    get_verts(g::AbstractSimpleGraph)

TBW
"""
get_verts(g::AbstractSimpleGraph) = deepcopy(g.verts)

function set_verts!(g::Graph, verts::Vector{<:Integer})
    g.verts = verts
end

"""
    get_edges(g::AbstractSimpleGraph)

TBW
"""
get_edges(g::AbstractSimpleGraph) = deepcopy(g.edges)

"""
    set_edges!(g::AbstractSimpleGraph, edges::Vector{<:Vector{<:Integer}})

TBW
"""
function set_edges!(g::AbstractSimpleGraph, edges::Vector{<:Vector{<:Integer}})
    if setdiff(vcat(edges...), get_verts(g)) != []
        error("Edges need to consist of vertices of the graph.")
    end
    g.edges = edges
end

"""
    no_concomponents(g::AbstractSimpleGraph)

TBW
"""
function no_concomponents(g::AbstractSimpleGraph)
    verts = collect(1:max(vcat(get_edges(g)...)...))
    sol = 0

    visited = [false for v in verts]

    function update!(k)
        visited[k] = true

        inc_edges = filter(e -> k in e, get_edges(g))
        if inc_edges == []
            return
        end

        adj = filter(l -> !visited[l], union(inc_edges...))
        for l in adj
            update!(l)
        end
    end

    for k in verts
        if !visited[k]
            update!(k)
            sol = sol + 1
        end
    end

    return sol
end


################## Embedded Graphs

abstract type AbstractEmbeddedGraph{S<:Real, T<:Integer} <:AbstractSimpleGraph{T} end

mutable struct Framework{S<:Real, T<:Integer} <:AbstractEmbeddedGraph{S,T}
    verts::Matrix{S}
    edges::Vector{Vector{T}}

    """
    Framework(verts::Matrix{<:Real}, edges::Vector{<:Vector{<:Integer}})

    TBW
    """
    function Framework(verts::AbstractMatrix{<:Real}, edges::Vector{<:Vector{<:Integer}})
        combgraph = Graph(edges=edges)
        if no_concomponents(combgraph) != 1
            error("Graph is not connected.")
        end
        if size(verts)[2] < max(vcat(edges...)...)
            error("Not every vertex of an edge has a coordinate assigned to it.")
        end

        S = typeof(verts[1,1])
        T = typeof(edges[1][1])
        return new{S,T}(verts, get_edges(combgraph))
    end
    
    function Framework(g::Graph)
        verts = rand(Float64, (3, max(get_verts(g)...)))
        return Framework(verts, get_edges(g))
    end
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
        """$(typeof(g)) in $(dimension(g))-space with $(length(get_verts(g))) vertices and $(length(get_edges(g))) edges.
        Edges: $(get_edges(g))
        """
    )
end

function edge_lengths(f::Framework)
    [sqrt(sqdist(get_verts(f)[:, e[1]], get_verts(f)[:, e[2]])) for e in get_edges(f)]
end