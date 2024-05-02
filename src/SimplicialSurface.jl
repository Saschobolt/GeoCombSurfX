import Graphs.connected_components, Graphs.SimpleGraph
import Base.Multimedia.display
using StaticArrays

include("Framework.jl")
include("Polyhedron.jl")
include("decomposition.jl")

################################################################################################################################################
###################################################### Type decalarations and constructors ######################################################
################################################################################################################################################

abstract type AbstractCombSimplicialSurface{T<:Integer} <: AbstractCombPolyhedron{T} end
abstract type AbstractSimplicialSurface{S<:Real,T<:Integer} <: AbstractPolyhedron{S,T} end

AbstractEmbOrCombSimplicialSurface{S<:Real,T<:Integer} = Union{AbstractSimplicialSurface{S,T},AbstractCombSimplicialSurface{T}}


mutable struct CombSimplicialSurface{T<:Integer} <: AbstractCombSimplicialSurface{T}
    verts::Vector{T}
    edges::Vector{MVector{2,T}}
    facets::Vector{MVector{3,T}}
    halfedges::Dict{MVector{2,T},HalfEdge{MVector{3,T}}}

    function CombSimplicialSurface(verts::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}})
        # @assert all([length(f) == 3 for f in facets]) "Facets of simplicial Surfaces are triangles."
        @assert all([length(filter(f -> issubset(e, f), facets)) in [1, 2] for e in edges]) "Edge of Simplicial Surface has to be edge of at least 1 and at most 2 facets."

        graph = Graphs.SimpleGraph(Graphs.Edge.(Tuple.(edges)))
        @assert length(connected_components(graph)) == 1 "1-Skeleton of Simplicial Surface has to be a connected graph."
        T = typeof(verts[1])
        surf = new{T}(verts, edges, facets)

        orient_facets!(surf)
        set_halfedges!(surf)

        return surf
    end
end

function CombSimplicialSurface(; verts=nothing, edges=nothing, facets=nothing)
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

mutable struct SimplicialSurface{S<:Real,T<:Integer} <: AbstractSimplicialSurface{S,T}
    verts::Matrix{S} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{MVector{2,T}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{MVector{3,T}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
    halfedges::Dict{MVector{2,T},HalfEdge{MVector{3,T}}}

    function SimplicialSurface(verts::AbstractMatrix{<:Real}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}}) # constructor
        comb_surf = CombSimplicialSurface(edges=edges, facets=facets)

        S = typeof(verts[1, 1])
        T = typeof(edges[1][1])
        return new{S,T}(verts, edges, get_facets(comb_surf), deepcopy(comb_surf.halfedges)) # in construction of comb_surf, the facets are oriented already
    end

    function SimplicialSurface(verts::AbstractVector{<:AbstractVector{<:Real}}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}})
        vertsmat = hcat(verts...)
        return SimplicialSurface(vertsmat, edges, facets)
    end

    function SimplicialSurface(; verts=nothing, edges=nothing, facets=nothing)
        comb_surf = CombSimplicialSurface(edges=edges, facets=facets)

        if isnothing(verts)
            n = length(get_verts(comb_surf))
            verts = rand(Float64, 3, n)
        end

        return SimplicialSurface(verts, get_edges(comb_surf), get_facets(comb_surf)) # in construction of comb_surf, the facets are oriented already
    end
end

function SimplicialSurface(poly::AbstractPolyhedron)
    poly_triang = triangulate(poly)
    return SimplicialSurface(get_verts(poly_triang), get_edges(poly_triang), get_facets(poly_triang))
end

function SimplicialSurface(comb_surf::AbstractCombSimplicialSurface)
    return SimplicialSurface(rand(Float64, 3, length(comb_surf.verts)), get_edges(comb_surf), get_facets(comb_surf))
end

function Polyhedron(surf::AbstractSimplicialSurface; atol::Real=1e-8, check_consistency::Bool=true)
    return Polyhedron(surf.verts, surf.edges, surf.facets; atol=atol, check_consistency=check_consistency)
end

function Framework(surf::AbstractSimplicialSurface)
    return Framework(surf.verts, surf.edges)
end

function CombSimplicialSurface(surf::AbstractSimplicialSurface)
    return CombSimplicialSurface(edges=get_edges(surf), facets=get_facets(surf))
end

################################################################################################################################################
###################################################### Simplicial Surface querying #############################################################
################################################################################################################################################
function vertex_degree(surf::AbstractEmbOrCombSimplicialSurface, v::Integer)
    return length(incfacets(surf, v))
end

characteristic(surf::AbstractCombSimplicialSurface) = length(surf.verts) - length(surf.edges) + length(surf.facets)
characteristic(surf::AbstractSimplicialSurface) = size(surf.verts)[2] - length(surf.edges) + length(surf.facets)

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

function iscactus(surf::AbstractSimplicialSurface)
    return iscactus(CombSimplicialSurface(surf))
end

"""
    cactus_distance_greedy(surf::AbstractCombSimplicialSurface)

TBW
"""
function cactus_distance_greedy(surf::AbstractCombSimplicialSurface)
    edgeturns = SVector{2,Int}[]

end

"""
    is_SimplicialSurface(poly::AbstractPolyhedron)

Determine wheter the Polyhedron poly is a Simplicial Surface.
"""
function is_SimplicialSurface(poly::AbstractEmbOrCombPolyhedron)
    try
        surf = CombSimplicialSurface(edges=get_edges(poly), facets=get_facets(poly))
    catch
        return false
    end
    return true
end

################################################################################################################################################
###################################################### Simplicial Surface manipulation and construction of new surfaces ########################
################################################################################################################################################

function set_facets!(surf::AbstractCombSimplicialSurface, facets::AbstractVector{<:AbstractVector{<:Integer}})
    if any([length(f) != 3 for f in facets])
        error("Facets of Simplicial Surfaces consist of 3 vertices")
    end

    # TODO: Meher asserts durchführen
    surf.facets = facets
end

function remove_vertex!(surf::AbstractCombSimplicialSurface, v::Integer; update_he::Bool=true)
    vertexmap = w -> w > v ? w - 1 : w
    setdiff!(surf.verts, [v])
    surf.verts = vertexmap.(surf.verts)
    setdiff!(surf.facets, incfacets(surf, v))
    surf.facets = [vertexmap.(f) for f in surf.facets]
    filter!(e -> !(v in e), surf.edges)
    surf.edges = [vertexmap.(e) for e in surf.edges]

    if update_he
        set_halfedges!(surf)
    end

    return surf
end

"""
    insert_butterfly!(surf::AbstractCombSimplicialSurface, edge1::AbstractVector{<:Integer}, edge2::AbstractVector{<:Integer}; is_oriented:Bool = false)


    Insert a butterfly along the vertex edge path described by edge1 and edge2. If is_oriented = true, the surface is assumed to be oriented.
"""
function insert_butterfly!(surf::AbstractCombSimplicialSurface, edge1::AbstractVector{<:Integer}, edge2::AbstractVector{<:Integer}; is_oriented::Bool=false)
    if !is_oriented
        orient_facets!(surf)
    end

    n = max(get_verts(surf)...)

    # vertex that is split by butterfly insertion
    split_vertex = Base.intersect(edge1, edge2)
    @assert length(split_vertex) == 1 "edge1 ($edge1) and edge2 ($edge2) have to have exactly one vertex in common, but intersection is $split_vertex"
    split_vertex = split_vertex[1]

    # edges should form a vertex edge path, so e1 = [v1, split_vertex], e2 = [split_vertex, v2]
    if indexin([split_vertex], edge1)[1] == 1
        e1 = reverse(edge1)
    else
        e1 = copy(edge1)
    end

    if indexin([split_vertex], edge2)[1] == 2
        e2 = reverse(edge2)
    else
        e2 = copy(edge2)
    end

    adj1 = adjfacets(surf, e1)
    @assert length(adj1) == 2 "edge1 ($edge1) has to have two adjacent facets, but has $(length(adj1)) ($adj1)."
    adj2 = adjfacets(surf, e2)
    @assert length(adj2) == 2 "edge2 ($edge2) has to have two adjacent facets, but has $(length(adj2)) ($adj2)."

    if length(Base.intersect(adj1, adj2)) == 1
        # first case: there are 3 adjacent facets in total -> wlog substitute split_vertex with n+1 in the adjacent facet, that contains both edges.
        f1 = Base.intersect(adj1, adj2)[1]

        # add new oriented facets to surf.
        if edge_direction(e1, f1) == 1
            push!(surf.facets, MVector{3}([split_vertex, n + 1, e1[1]]))
            push!(surf.facets, MVector{3}([split_vertex, e2[2], n + 1]))
        else
            push!(surf.facets, MVector{3}([split_vertex, e1[1], n + 1]))
            push!(surf.facets, MVector{3}([split_vertex, n + 1, e2[2]]))
        end

        surf.facets[indexin([f1], surf.facets)[1]][indexin([split_vertex], f1)[1]] = n + 1
    else
        # second case: there are 4 adjacent facets in total -> wlog substitute split_vertex with n+1 in the adjacent facet of e1, where e1 is facing forwards, and in the adjacent facet of e2, where e2 is facing forwards.
        f1 = filter(f -> edge_direction(e1, f) == 1, adj1)[1]
        f2 = filter(f -> edge_direction(e2, f) == 1, adj2)[1]

        # get gallery connecting f1 to f2, which has to be updated.
        gal_facets = [f1]
        gal_edges = filter(e -> !(e in [e1, reverse(e1)]) && split_vertex in e, incedges(surf, f1))
        while gal_facets[end] != f2
            next = filter(f -> (split_vertex in f) && !issubset(e1, f), setdiff(adjfacets(surf, gal_facets[end]), gal_facets))[1]
            append!(gal_edges, filter(e -> !(e in [e2, reverse(e2)]) && split_vertex in e, setdiff(incedges(surf, next), gal_edges)))
            push!(gal_facets, next)
        end

        # set split_vertex to n+1 in all facets in gal_facets and all edges in gal_edges
        setdiff!(surf.facets, gal_facets)
        setdiff!(surf.edges, gal_edges)

        map!(f -> setindex!(f, n + 1, indexin([split_vertex], f)[1]), gal_facets, gal_facets)
        map!(e -> setindex!(e, n + 1, indexin([split_vertex], e)[1]), gal_edges, gal_edges)

        append!(surf.facets, gal_facets)
        append!(surf.edges, gal_edges)

        # surf.facets[indexin([f1], surf.facets)[1]][indexin([split_vertex], f1)[1]] = n+1
        # surf.facets[indexin([f2], surf.facets)[1]][indexin([split_vertex], f2)[1]] = n+1

        # add new oriented facets to surf.
        push!(surf.facets, MVector{3}([split_vertex, n + 1, e1[1]]))
        push!(surf.facets, MVector{3}([split_vertex, e2[2], n + 1]))
    end

    # add new edges to surf.
    push!(surf.edges, MVector{2}([split_vertex, n + 1]))
    push!(surf.edges, MVector{2}([e1[1], n + 1]))
    push!(surf.edges, MVector{2}([n + 1, e2[2]]))

    # add new vertex to surf
    push!(surf.verts, n + 1)

    # update halfedges
    set_halfedges!(surf, is_oriented=true)

    return surf
end

function insert_butterfly!(surf::AbstractSimplicialSurface, edge1::AbstractVector{<:Integer}, edge2::AbstractVector{<:Integer}, p::Union{AbstractVector,Nothing}=nothing; is_oriented::Bool=false)
    comb_surf = CombSimplicialSurface(edges=get_edges(surf), facets=get_facets(surf))
    insert_butterfly!(comb_surf, edge1, edge2, is_oriented=is_oriented)

    if isnothing(p)
        surf.verts = hcat(surf.verts, rand(Float64, 3, 1))
    else
        surf.verts = hcat(surf.verts, p)
    end
    surf.edges = get_edges(comb_surf)
    surf.facets = get_facets(comb_surf)

    return surf
end

"""
    insert_butterfly!(surf::AbstractEmbOrCombSimplicialSurface, edge1::AbstractVector{<:Integer}, edge2::AbstractVector{<:Integer}, p::Union{AbstractVector, Nothing} = nothing; is_oriented::Bool = false)

Insert a butterfly along the vertex edge path described by edge1 and edge2. If is_oriented = true, the surface is assumed to be oriented. 
If surf is an embedded simplicial surface, the new vertex is placed at p. If p is nothing, a random point is chosen. 
If surf is a combinatorial simplicial surface, p has to be nothing.
"""
function insert_butterfly!(surf::AbstractEmbOrCombSimplicialSurface, edge1::AbstractVector{<:Integer}, edge2::AbstractVector{<:Integer}, p::Union{AbstractVector,Nothing}=nothing; is_oriented::Bool=false)
    if isnothing(p)
        return insert_butterfly!(surf, edge1, edge2, is_oriented=is_oriented)
    else
        @assert typeof(surf) <: AbstractSimplicialSurface "p has to be nothing if surf is a combinatorial simplicial surface, but got $p."
        return insert_butterfly!(surf, edge1, edge2, p, is_oriented=is_oriented)
    end
end

"""
    insert_butterfly(surf::AbstractCombSimplicialSurface, edge1::AbstractVector{<:Integer}, edge2::AbstractVector{<:Integer}; is_oriented:Bool = false)

Insert a butterfly along the vertex edge path described by edge1 and edge2. If is_oriented = true, the surface is assumed to be oriented.
"""
function insert_butterfly(surf::AbstractCombSimplicialSurface, edge1::AbstractVector{<:Integer}, edge2::AbstractVector{<:Integer}; is_oriented::Bool=false)
    surf_copy = deepcopy(surf)
    insert_butterfly!(surf_copy, edge1, edge2, is_oriented=is_oriented, update_he=true)

    return surf_copy
end

"""
    random_simplsphere(n::Integer)

Construct a random simplicial sphere with n vertices. The sphere is constructed by starting with a tetrahedron and applying random butterfly insertions.
"""
function random_simplsphere(n::Integer)
    sphere = CombSimplicialSurface(verts=[1, 2, 3, 4], edges=[[1, 2], [2, 3], [3, 1], [1, 4], [2, 4], [3, 4]], facets=[[1, 2, 3], [4, 2, 1], [4, 3, 2], [1, 3, 4]])
    for _ in 1:n-4
        v = rand(1:length(sphere.verts))
        e1 = rand(incedges(sphere, v))
        e2 = rand(setdiff(incedges(sphere, v), [e1]))
        insert_butterfly!(sphere, e1, e2; is_oriented=true)
    end

    return sphere
end

random_emb_simplsphere(n::Integer) = SimplicialSurface(random_simplsphere(n))

# function isadjacent(surf::AbstractEmbOrCombSimplicialSurface, facetoredge::AbstractVector{<:Integer}, facet::AbstractVector{<:Integer}; check::Bool=true)
#     if check
#         @assert facetoredge in union(get_edges(surf), get_facets(surf)) || reverse(facetoredge) in union(get_edges(surf), get_facets(surf)) "facetoredge has to be an edge or facet of the surface, but got $(facetoredge)."
#         @assert facet in get_facets(surf) || reverse(facet) in get_facets(surf) "facet has to be a facet of the surface, but got $(facet)."
#     end

#     return length(Base.intersect(facetoredge, facet)) == 2
# end

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

    # update halfedges
    set_halfedges!(surf, is_oriented=false)

    return surf
end

function edgeturn(surf::AbstractEmbOrCombSimplicialSurface, e::AbstractVector{<:Integer})
    surf_copy = deepcopy(surf)
    edgeturn!(surf_copy, e)

    return surf_copy
end

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

    set_halfedges!(surf)

    return surf
end

function append_tetrahedron!(surf::AbstractCombSimplicialSurface, f::AbstractVector{<:Integer}; check::Bool=true, update_he::Bool=true)
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
    append!(surf.facets, [[facet[1], facet[2], n + 1], [facet[2], facet[3], n + 1], [facet[3], facet[1], n + 1]])
    setdiff!(surf.facets, [facet])
    append!(surf.edges, [MVector{2}([x, n + 1]) for x in facet])
    push!(surf.verts, n + 1)

    if update_he
        set_halfedges!(surf, is_oriented=true)
    end

    return surf
end

function append_tetrahedron!(surf::AbstractSimplicialSurface, f::AbstractVector{<:Integer}, p::Union{Nothing,AbstractVector{<:Real}}=nothing; check::Bool=true, update_he::Bool=true)
    comb_surf = CombSimplicialSurface(edges=get_edges(surf), facets=get_facets(surf))
    append_tetrahedron!(comb_surf, f, check=check, update_he=update_he)

    if isnothing(p)
        surf.verts = hcat(surf.verts, rand(Float64, 3, 1))
    else
        surf.verts = hcat(surf.verts, p)
    end

    surf.edges = get_edges(comb_surf)
    surf.facets = get_facets(comb_surf)

    return surf
end

"""
    append_tetrahedron!(surf::AbstractEmbOrCombSimplicialSurface, f::AbstractVector{<:Integer}, p::Union{Nothing, AbstractVector{<:Real}} = nothing; check::Bool = true)

Append a tetrahedron to surf by removing facet f and adding a tetrahedron in its place. 
In case surf is an embedded surface, the new vertex is placed at p. If p is nothing, a random point is chosen. If surf is a combinatorial simplicial surface, p has to be nothing.
If check is set to true, it is checked that f is a facet of surf.
"""
function append_tetrahedron!(surf::AbstractEmbOrCombSimplicialSurface, f::AbstractVector{<:Integer}, p::Union{Nothing,AbstractVector{<:Real}}=nothing; check::Bool=true, update_he::Bool=true)
    if isnothing(p)
        return append_tetrahedron!(surf, f, check=check)
    else
        @assert typeof(surf) <: AbstractSimplicialSurface "p has to be nothing if surf is a combinatorial simplicial surface, but got $p."
        return append_tetrahedron!(surf, f, p, check=check, update_he=update_he)
    end
end

function random_cactus(n::Integer; colored::Bool=false)
    cactus = CombSimplicialSurface(verts=[1, 2, 3, 4], edges=[[1, 2], [2, 3], [3, 1], [1, 4], [2, 4], [3, 4]], facets=[[1, 2, 3], [4, 2, 1], [4, 3, 2], [1, 3, 4]])
    if colored
        color_dict = Dict([1, 2] => 1, [2, 3] => 2, [3, 1] => 3, [1, 4] => 2, [2, 4] => 3, [3, 4] => 1)
        cactus = ColoredSimplicialSurface(cactus, color_dict)
    end

    for _ in 1:n-1
        i = rand(1:length(cactus.facets))
        append_tetrahedron!(cactus, cactus.facets[i], check=false)
    end

    return cactus
end

function random_emb_cactus(n::Integer)
    cactus = random_cactus(n)

    return SimplicialSurface(cactus)
end

################################################################################################################################################
###################################################### Coloured Simplicial Surfaces ############################################################
################################################################################################################################################
abstract type AbstractColoredSimplicialSurface{T<:Integer} <: AbstractCombSimplicialSurface{T} end

mutable struct ColoredSimplicialSurface{T<:Integer} <: AbstractColoredSimplicialSurface{T}
    verts::Vector{T}
    edges::Vector{MVector{2,T}}
    facets::Vector{MVector{3,T}}
    halfedges::Dict{MVector{2,T},HalfEdge{MVector{3,T}}}
    color_dict::Dict{MVector{2,T},Int} # Dictionary containing edge colors

    function ColoredSimplicialSurface(surf::AbstractCombSimplicialSurface{T}, color_dict::Dict{<:AbstractVector{<:Integer},<:Integer}) where {T<:Integer}
        @assert issubset(surf.edges, keys(color_dict)) "color_dict has to have edges of surf as keys."

        return new{T}(surf.verts, surf.edges, surf.facets, surf.halfedges, color_dict)
    end
end

function ColoredSimplicialSurface(verts, edges, facets, color_dict)
    surf = CombSimplicialSurface(verts, edges, facets)
    return ColoredSimplicialSurface(surf, color_dict)
end

function ColoredSimplicialSurface(; verts, edges, facets, color_dict)
    surf = CombSimplicialSurface(verts=verts, edges=edges, facets=facets)
    return ColoredSimplicialSurface(surf, color_dict)
end

function display(surf::AbstractColoredSimplicialSurface)
    print("""$(typeof(surf)) with $(length(get_verts(surf))) vertices, $(length(get_edges(surf))) edges, $(length(get_facets(surf))) facets, $(length(unique(collect(values(surf.color_dict))))) colors.\n 
    Edges:  $(get_edges(surf))
    Facets: $(get_facets(surf))
    Colors: $(colors(surf)) \n""")
end

function colors(surf::AbstractColoredSimplicialSurface)
    return sort(unique(collect(values(surf.color_dict))))
end

function color(surf::AbstractColoredSimplicialSurface, e::AbstractVector{<:Integer})
    if e in keys(surf.color_dict)
        return surf.color_dict[e]
    elseif reverse(e) in keys(surf.color_dict)
        return surf.color_dict[reverse(e)]
    end

    error("Either e ($e) is not an edge of the surface surf, or it is not colored.")
end

function colortype(surf::AbstractColoredSimplicialSurface, f::AbstractVector{<:Integer}; check::Bool=true)
    if check
        @assert Set(f) in Set.(surf.facets)
    end

    return [color(surf, MVector{2}([f[1], f[2]])), color(surf, MVector{2}([f[2], f[3]])), color(surf, MVector{2}([f[3], f[1]]))]
end

function congruencetypes(surf::AbstractColoredSimplicialSurface)
    return unique(sort.([colortype(surf, f) for f in surf.facets]))
end

"""
    is_congcolored(surf::AbstractColoredSimplicialSurface)

Return whether the colored simplicial surface surf corresponds to a simplicial surface with congruent triangles, i.e. every facet has the same colors on the edges.
"""
function is_congcolored(surf::AbstractColoredSimplicialSurface)
    return length(congruencetypes(surf)) == 1
end

function append_tetrahedron!(surf::AbstractColoredSimplicialSurface, f::AbstractVector{<:Integer}; check::Bool=true, update_he::Bool=true, color_new_edges::Bool=true)
    invoke(append_tetrahedron!, Tuple{AbstractCombSimplicialSurface,AbstractVector{<:Integer}}, surf, f, check=check, update_he=update_he)


    if color_new_edges
        for i in 1:3
            e1 = f[[i, mod1(i + 1, 3)]]
            e2 = f[[i, mod1(i + 2, 3)]]
            col = setdiff(colors(surf), [color(surf, e1), color(surf, e2)])[1]
            surf.color_dict[[f[i], surf.verts[end]]] = col
        end
    end
    return surf
end

"""
    edge_type(surf::AbstractColoredSimplicialSurface, e::AbstractVector{<:Integer})

Return the edge type of an edge of a colored simplicial surface. Returns "r" for a rotation edge, "m" for a mirror edge, "mr" for a simultaneous mirror and reflection edge. 
Returns "b" if edge is a boundary edge and 0, if edge is an interior edge but neither "r" nor "m".
"""
function edge_type(surf::AbstractColoredSimplicialSurface, e::AbstractVector{<:Integer})
    he = surf.halfedges[e]

    # if e is a boundary edge, return 0
    if isnothing(he.twin)
        return "b"
    end

    # if adjacent facets of e are not congruent, return 0
    if sort(colortype(surf, he.face)) != sort(colortype(surf, he.twin.face))
        return "0"
    end

    # if adjacent facets of e are equilateral, return "mr", as 
    if length(unique(colortype(surf, he.face))) == 1
        return "mr"
    end

    e1 = [e[1], setdiff(he.face, e)[1]]
    e2 = [e[1], setdiff(he.twin.face, e)[1]]

    if color(surf, e1) == color(surf, e2)
        return "m"
    end

    return "r"
end

function is_tamecolored(surf::AbstractColoredSimplicialSurface)
    edge_types = Dict{Int,String}()

    for e in surf.edges
        type = edge_type(surf, e)
        if type == "0"
            return false
        end

        if type == "b"
            continue
        end

        col = color(surf, e)
        ref = get!(edge_types, col, type)
        if ref != type
            return false
        end
    end

    return true
end