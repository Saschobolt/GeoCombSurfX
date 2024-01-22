include("Polyhedron.jl")
include("decomposition.jl")
############################################################################
#       Combinatorial Simplicial Surfaces
############################################################################
abstract type AbstractCombSimplicialSurface{T<:Integer} <: AbstractSimpleGraph{T} end

mutable struct CombSimplicialSurface{T<:Integer} <: AbstractCombSimplicialSurface{T}
    verts::Vector{T}
    edges::Vector{Vector{T}}
    facets::Vector{Vector{T}}

    function CombSimplicialSurface(verts::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}})
        @assert all([length(f) == 3 for f in facets]) "Facets of simplicial Surfaces are triangles."
        @assert all([length(filter(f -> issubset(e, f), facets)) in [1,2] for e in edges]) "Edge of Simplicial Surface has to be edge of at least 1 and at most 2 facets."

        graph = Graph(verts, edges)
        @assert no_concomponents(graph) == 1 "1-Skeleton of Simplicial Surface has to be a connected graph."
        T = typeof(verts[1])

        poly = Polyhedron(rand(Float64, 3, length(verts)), edges, facets)
        orient_facets!(poly)
        return new{T}(verts, edges, get_facets(poly))
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

        graph = Graph(verts = verts, edges = edges)
        return CombSimplicialSurface(get_verts(graph), get_edges(graph), facets)
    end
end


get_facets(surf::AbstractCombSimplicialSurface) = deepcopy(surf.facets)

function set_facets!(surf::AbstractCombSimplicialSurface, facets::AbstractVector{<:AbstractVector{<:Integer}})
    if any([length(f) != 3 for f in facets])
        error("Facets of Simplicial Surfaces consist of 3 vertices")
    end

    # TODO: Meher asserts durchführen
    surf.facets = facets
end

function orient_facets!(surf::AbstractCombSimplicialSurface)
    new_surf = CombSimplicialSurface(get_verts(surf), get_edges(surf), get_facets(surf)) # in construction of CombSimplicialSurface the facets are oriented
    set_facets!(surf, get_facets(new_surf))
end

function orient_facets(surf::AbstractCombSimplicialSurface)
    surf_copy = deepcopy(surf)
    
    orient_facets!(surf_copy)
    return surf_copy
end

function boundary(surf::AbstractCombSimplicialSurface)
    return filter(e -> length(filter(f -> issubset(e, f), get_facets(surf))) == 1, get_edges(surf))
end

function incedges(surf::AbstractCombSimplicialSurface, f::AbstractVector{<:Integer})
    return incedges(SimplicialSurface(surf), f)
end

function incedges(surf::AbstractCombSimplicialSurface, v::Int)
    return incedges(SimplicialSurface(surf), v)
end

function incfacets(surf::AbstractCombSimplicialSurface, vertexarray::AbstractVector{<:Integer})
    return incfacets(SimplicialSurface(surf), vertexarray)    
end

function incfacets(surf::AbstractCombSimplicialSurface, v::Integer)
    return incfacets(SimplicialSurface(surf), v)    
end

function edgeturn!(surf::AbstractCombSimplicialSurface, e::AbstractVector{<:Integer})
    if !(e in get_edges(surf) || reverse(e) in get_edges(surf))
        error("e has to be an edge of surf, but got $(e).")
    end
    butterfly = incfacets(surf, e)
    tips = symdiff(butterfly...)

    # delete e from edges and add tips as new edge
    edges = get_edges(surf)
    if tips in edges
        error("The tips of the butterfly $(butterfly) having $(e) as an interior edge is itself an edge of the Surface ($(tips)). Thus the edge is not turnable.")
    end
    edges = setdiff(edges, [e, reverse(e)])
    push!(edges, tips)

    set_edges!(surf, edges)

    # delete old facets and append new ones
    facets = get_facets(surf)
    facets = setdiff(facets, butterfly)

    newfacets = [[tips[1], tips[2], e[1]], [tips[1], tips[2], e[2]]]
    append!(facets, newfacets)
    set_facets!(surf, facets)

    orient_facets!(surf)
end

function edgeturn(surf::AbstractCombSimplicialSurface, e::AbstractVector{<:Integer})
    surf_copy = deepcopy(surf)
    edgeturn!(surf_copy, e)

    return surf_copy
end



############################################################################
#       Embedded Simplicial Surfaces
############################################################################
abstract type AbstractSimplicialSurface{S<:Real, T<:Integer} <: AbstractPolyhedron{S,T} end

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

function CombSimplicialSurface(surf::AbstractSimplicialSurface)
    return CombSimplicialSurface(edges = get_edges(surf), facets = get_facets(surf))
end

"""
    is_SimplicialSurface(poly::AbstractPolyhedron)

Determine wheter the Polyhedron poly is a Simplicial Surface.
"""
function is_SimplicialSurface(poly::AbstractPolyhedron)
    try
        surf = CombSimplicialSurface(edges = get_edges(poly), facets = get_facets(poly))
    catch
        return false
    end
    return true
end

function boundary(surf::AbstractSimplicialSurface)
    return boundary(CombSimplicialSurface(surf))
end

function edgeturn!(surf::AbstractSimplicialSurface, e::AbstractVector{<:Integer})
    comb_surf = CombSimplicialSurface(surf)
    edgeturn!(comb_surf, e)

    set_edges!(surf, get_edges(comb_surf))
    set_facets!(surf, get_facets(comb_surf))
end

function edgeturn(surf::AbstractSimplicialSurface, e::AbstractVector{<:Integer})
    surf_copy = deepcopy(surf)
    edgeturn!(surf_copy, e)
    return surf_copy
end