import Graphs.connected_components, Graphs.SimpleGraph
using StaticArrays

include("Framework.jl")
include("Polyhedron.jl")
include("decomposition.jl")

abstract type AbstractCombSimplicialSurface{T<:Integer} <: AbstractCombPolyhedron{T} end
abstract type AbstractSimplicialSurface{S<:Real, T<:Integer} <: AbstractPolyhedron{S,T} end

AbstractEmbOrCombSimplicialSurface{S<:Real, T<:Integer} = Union{AbstractSimplicialSurface{S,T}, AbstractCombSimplicialSurface{T}}


mutable struct CombSimplicialSurface{T<:Integer} <: AbstractCombSimplicialSurface{T}
    verts::Vector{T}
    edges::Vector{SVector{2, T}}
    facets::Vector{SVector{3, T}}

    function CombSimplicialSurface(verts::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}})
        # @assert all([length(f) == 3 for f in facets]) "Facets of simplicial Surfaces are triangles."
        @assert all([length(filter(f -> issubset(e, f), facets)) in [1,2] for e in edges]) "Edge of Simplicial Surface has to be edge of at least 1 and at most 2 facets."

        graph = Graphs.SimpleGraph(Graphs.Edge.(Tuple.(edges)))
        @assert length(connected_components(graph)) == 1 "1-Skeleton of Simplicial Surface has to be a connected graph."
        T = typeof(verts[1])

        poly = Polyhedron(rand(Float64, 3, length(verts)), edges, facets)
        orient_facets!(poly)
        return new{T}(verts, SVector{2}.(edges), SVector{3}.(get_facets(poly)))
    end

    function CombSimplicialSurface(; verts = nothing, edges = nothing, facets = nothing)
        if isnothing(edges) && !isnothing(facets)
            edges = Vector{Int64}[]
            for f in facets
                e1 = [f[1], f[2]]
                e2 = [f[2], f[3]]
                e3 = [f[3], f[1]]

                append!(edges, [e1, e2, e3])
            end
            edges = collect.(unique(Set.(edges)))
        elseif isnothing(facets) && !isnothing(edges)
            error("not implemented yet.")
        end

        if isnothing(verts)
            n = max(vcat(Vector.(edges)...)...)
            verts = collect(1:n)
        end

        return CombSimplicialSurface(verts, edges, facets)
    end
end

function set_facets!(surf::AbstractCombSimplicialSurface, facets::AbstractVector{<:AbstractVector{<:Integer}})
    if any([length(f) != 3 for f in facets])
        error("Facets of Simplicial Surfaces consist of 3 vertices")
    end

    # TODO: Meher asserts durchführen
    surf.facets = facets
end

function remove_vertex!(surf::AbstractCombSimplicialSurface, v::Integer)
    vertexmap = w -> w > v ? w - 1 : w
    setdiff!(surf.verts, [v])
    surf.verts = vertexmap.(surf.verts)
    setdiff!(surf.facets, incfacets(surf, v))
    surf.facets = [vertexmap.(f) for f in surf.facets]
    filter!(e -> !(v in e), surf.edges)
    surf.edges = [vertexmap.(e) for e in surf.edges]

    return surf
end

mutable struct SimplicialSurface{S<:Real, T<:Integer} <: AbstractSimplicialSurface{S,T}
    verts::Matrix{S} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{Vector{T}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{Vector{T}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
    function SimplicialSurface(verts::AbstractMatrix{<:Real}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}}) # constructor
        comb_surf = CombSimplicialSurface(edges = edges, facets = facets)

        S = typeof(verts[1,1])
        T = typeof(edges[1][1])
        return new{S,T}(verts, edges, get_facets(comb_surf)) # in construction of comb_surf, the facets are oriented already
    end

    function SimplicialSurface(verts::AbstractVector{<:AbstractVector{<:Real}}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}})
        vertsmat = hcat(verts...)
        return SimplicialSurface(vertsmat, edges, facets)
    end
    
    function SimplicialSurface(;verts = nothing, edges = nothing, facets = nothing)    
        comb_surf = CombSimplicialSurface(edges = edges, facets = facets)

        if isnothing(verts)
            n = length(get_verts(comb_surf))
            verts = rand(Float64, 3, n)
        end
        S = typeof(verts[1,1])
        T = typeof(get_edges(comb_surf)[1][1])
        return SimplicialSurface(verts, get_edges(comb_surf), get_facets(comb_surf)) # in construction of comb_surf, the facets are oriented already
    end
end

function SimplicialSurface(poly::AbstractPolyhedron)
    poly_triang = triangulate(poly)
    return SimplicialSurface(get_verts(poly_triang), get_edges(poly_triang), get_facets(poly_triang))
end

function SimplicialSurface(comb_surf::AbstractCombSimplicialSurface)
    return SimplicialSurface(edges = get_edges(comb_surf), facets = get_facets(comb_surf))
end

function Polyhedron(surf::AbstractSimplicialSurface)
    return Polyhedron(surf.verts, surf.edges, surf.facets)
end

function Framework(surf::AbstractSimplicialSurface)
    return Framework(surf.verts, surf.edges)
end

function CombSimplicialSurface(surf::AbstractSimplicialSurface)
    return CombSimplicialSurface(edges = get_edges(surf), facets = get_facets(surf))
end

function remove_vertex!(surf::AbstractSimplicialSurface, v::Integer)
    vertexmap = w -> w > v ? w - 1 : w
    surf.verts = surf.verts[:, setdiff(1:end, v)]
    setdiff!(surf.facets, incfacets(surf, v))
    surf.facets = [vertexmap.(f) for f in surf.facets]
    filter!(e -> !(v in e), surf.edges)
    surf.edges = [vertexmap.(e) for e in surf.edges]

    return surf
end

function isadjacent(surf::AbstractEmbOrCombSimplicialSurface, facetoredge::AbstractVector{<:Integer}, facet::AbstractVector{<:Integer}; check::Bool = true)
    if check
        @assert facetoredge in union(get_edges(surf), get_facets(surf)) || reverse(facetoredge) in union(get_edges(surf), get_facets(surf)) "facetoredge has to be an edge or facet of the surface, but got $(facetoredge)."
        @assert facet in get_facets(surf) || reverse(facet) in get_facets(surf) "facet has to be a facet of the surface, but got $(facet)."
    end

    return length(Base.intersect(facetoredge, facet)) == 2
end

function edgeturn!(surf::AbstractEmbOrCombSimplicialSurface, e::AbstractVector{<:Integer})
    if !(e in get_edges(surf) || reverse(e) in get_edges(surf))
        error("e has to be an edge of surf, but got $(e).")
    end
    butterfly = incfacets(surf, e)
    tips = symdiff(butterfly...)

    # delete e from edges and add tips as new edge
    if tips in surf.edges || reverse(tips) in surf.edges
        error("The tips of the butterfly $(butterfly) having $(e) as an interior edge is itself an edge of the Surface ($(tips)). Thus the edge is not turnable.")
    end
    setdiff!(surf.edges, [e, reverse(e)])
    push!(surf.edges, tips)
    
    # delete old facets and append new ones
    setdiff!(surf.facets, butterfly)

    append!(surf.facets, [[e[1], tips[1], tips[2]], [e[2], tips[2], tips[1]]])

    return surf
end

function edgeturn(surf::AbstractEmbOrCombSimplicialSurface, e::AbstractVector{<:Integer})
    surf_copy = deepcopy(surf)
    edgeturn!(surf_copy, e)

    return surf_copy
end

function vertex_degree(surf::AbstractEmbOrCombSimplicialSurface, v::Integer)
    return length(incfacets(surf, v))
end

characteristic(surf::AbstractCombSimplicialSurface) = length(surf.verts) - length(surf.edges) + length(surf.facets)
characteristic(surf::AbstractSimplicialSurface) = size(surf.verts)[2] - length(surf.edges) + length(surf.facets)

"""
    remove_tetrahedron(surf::AbstractCombSimplicialSurface, v::Integer)

Remove vertex v of vertex degree 3 from the simplicial surface. Add the base of the correspnding tetrahedron as a facet of the surface.
"""
function remove_tetrahedron!(surf::AbstractEmbOrCombSimplicialSurface, v::Integer)
    tetrahedron = incfacets(surf, v)
    base = unique(setdiff(vcat(tetrahedron...), [v]))
    if length(base) > 3
        error("v ($v) is expected to be of vertex degree 3, but got incident facets $tetrahedron")
    end
    push!(surf.facets, base)
    remove_vertex!(surf, v)

    return surf
end


"""
    append_tetrahedron!(surf::AbstractCombSimplicialSurface, f::AbstractVector{<:Integer})

Append a tetrahedron to surf by removing facet f and adding a tetrahedron in its place. If check is set to true, it is checked that f is a facet of surf.
"""
function append_tetrahedron!(surf::AbstractCombSimplicialSurface, f::AbstractVector{<:Integer}; check::Bool = true)
    if check
        i = findfirst(x -> Base.intersect(x, f) == x, surf.facets)
        if isnothing(i)
            error("f needs to be a facet of surf, but got $(f)")
        end
        facet = deepcopy(surf.facets[i])
    else
        facet = f
    end

    n = length(surf.verts)
    append!(surf.facets, [push(e, n+1) for e in incedges(surf, facet)])
    setdiff!(surf.facets, [facet])
    append!(surf.edges, [SVector{2}([x, n+1]) for x in facet])
    push!(surf.verts, n+1)

    return surf
end

function random_cactus(n::Integer)
    cactus = CombSimplicialSurface(facets = [[1,2,3], [1,2,4], [2,3,4], [1,3,4]])
    for _ in 1:n-1
        i = rand(1:length(cactus.facets))
        append_tetrahedron!(cactus, cactus.facets[i], check = false)
    end

    return cactus
end

function iscactus(surf::AbstractCombSimplicialSurface)
    # Cacti have genus 0
    if characteristic(surf) != 2
        return false
    end

    cactus = deepcopy(surf)
    while length(cactus.verts) > 4
        v = get_verts(cactus)[findfirst(x -> vertex_degree(cactus, x) == 3, cactus.verts)]
        if isnothing(v)
            return false
        end
        remove_tetrahedron!(cactus, v)
    end

    return true
end

"""
    cactus_distance_greedy(surf::AbstractCombSimplicialSurface)

TBW
"""
function cactus_distance_greedy(surf::AbstractCombSimplicialSurface)
    edgeturns = SVector{2, Int}[]
    
end

"""
    is_SimplicialSurface(poly::AbstractPolyhedron)

Determine wheter the Polyhedron poly is a Simplicial Surface.
"""
function is_SimplicialSurface(poly::AbstractEmbOrCombPolyhedron)
    try
        surf = CombSimplicialSurface(edges = get_edges(poly), facets = get_facets(poly))
    catch
        return false
    end
    return true
end
