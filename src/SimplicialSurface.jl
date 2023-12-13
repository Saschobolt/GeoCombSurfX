include("Polyhedron.jl")
include("decomposition.jl")

abstract type AbstractSimplicialSurface{S<:Real, T<:Integer} <: AbstractPolyhedron{S,T} end

mutable struct SimplicialSurface{S<:Real, T<:Integer} <: AbstractSimplicialSurface{S,T}
    verts::Vector{Vector{S}} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{Vector{T}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{Vector{T}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
    function SimplicialSurface(verts::Vector{Vector{S}}, edges::Vector{Vector{T}}, facets::Vector{Vector{T}}) where S<:Real where T<:Integer # constructor
        @assert all([length(verts[i]) == length(verts[1]) for i in 1:length(verts)]) "Dimension is not well defined."
        @assert all([length(f) == 3 for f in facets]) "Facets of a simplicial surface are triangular and thus contain exactly three vertices."
        @assert all([1 <= count(f -> issubset(e,f), facets) <= 2 for e in edges]) "Every edge needs to be contained in at least one and at most two facets."
        # TODO: Umbrella condition als assert hinzufügen.

        return orient_facets(new{S,T}(verts, edges, facets))
    end
end

function SimplicialSurface(;verts::Vector{<:Vector{<:Real}} = nothing, edges::Vector{<:Vector{<:Integer}} = nothing, facets::Vector{<:Vector{<:Integer}})
    if isnothing(verts)

    end

    if isnothing(edges) && !isnothing(facets)

    elseif isnothing(facets) && !isnothing(edges)
        error("not implemented yet.")
    else

    edges = collect.(union([[Set([f[1], f[2]]), Set([f[2], f[3]]), Set([f[1], f[3]])] for f in facets]...)) # 2-element subsets of facets.
    return SimplicialSurface(verts, edges, facets)
end

function SimplicialSurface(verts::AbstractMatrix{<:Real}, edges::Vector{<:Vector{<:Integer}}, facets::Vector{<:Vector{<:Integer}})
    vertsarray = [verts[:,i] for i in 1:size(verts)[2]]
    return SimplicialSurface(vertsarray, edges, facets)
end

function SimplicialSurface(verts::AbstractMatrix{<:Real}, facets::Vector{<:Vector{<:Integer}})
    vertsarray = [verts[:,i] for i in 1:size(verts)[2]]
    return SimplicialSurface(vertsarray, facets)
end

function SimplicialSurface(poly::AbstractPolyhedron)
    polytriang = triangulate(poly)
    return SimplicialSurface(get_verts(polytriang), get_edges(polytriang), get_facets(polytriang))
end