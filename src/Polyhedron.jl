using GenericLinearAlgebra
using PolygonOps
import Base.Multimedia.display
import Base.==
import Graphs.connected_components

include("Framework.jl")
include("affine_geometry.jl")
include("polygonal_geometry.jl")

abstract type AbstractPolyhedron{S<:Real,T<:Integer} <: AbstractEmbeddedGraph{S, T} end
abstract type AbstractCombPolyhedron{T<:Integer} end

AbstractEmbOrCombPolyhedron = Union{AbstractPolyhedron, AbstractCombPolyhedron}

mutable struct Polyhedron{S<:Real, T<:Integer} <:AbstractPolyhedron{S, T}
    verts::Matrix{S} # matrix of coordinates of the vertices. Columns correspond to vertices.
    edges::Vector{Vector{T}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{Vector{T}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
    # TODO: Neuen constructor mit optionalen Argumenten (wenn nur coordinates gegeben werden, ist Ergebnis die konvexe Hülle der Punkte + check der Dimension
    # wenn nur Facets gegeben sind, werden Edges automatisch gesetzt und es wird gecheckt, dass Vertizes auf einer Facet koplanar aber nicht kollinear sind)
    function Polyhedron(verts::AbstractMatrix{<:Real}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}}; atol::Real = 1e-8, check_consistency::Bool = true)
        for f in facets
            if length(f) < 3
                error("Facets need to consist of at least 3 vertices.")
            elseif length(f) == 3
                continue
            elseif affinedim(verts[:,f]; atol = atol) != 2
                error("Facets with more than 3 vertices need to span a space of affine dimension 2.")
            end                            
        end

        if check_consistency
            for f in facets
                n = length(f)
                for i in 1:n
                    if !(f[[mod1(i, n), mod1(i+1, n)]] in edges || f[[mod1(i+1, n), mod1(i, n)]] in edges)
                        error("Facets and edges are not consistent.")
                    end
                end
            end
        end
        
        S = typeof(verts[1,1])
        T = typeof(edges[1][1])
        poly = new{S,T}(verts, edges, facets)

        if length(connected_components(SimpleGraph(poly))) > 1
            error("Skeleton of polyhedron is not connected.")
        end

        if check_consistency
            for e in edges
                if length(adjfacets(poly, e)) > 2
                    error("Edge $(e) has to be edge of at least one and at most two facets, but is adjacent to $(adjfacets(poly, e)).")
                end
            end
        end
        return orient_facets(poly)
    end

    function Polyhedron(verts::AbstractVector{<:AbstractVector{<:Real}}, edges::AbstractVector{<:AbstractVector{<:Integer}}, facets::AbstractVector{<:AbstractVector{<:Integer}}; atol::Real = 1e-8)
        return Polyhedron(hcat(verts...), edges, facets; atol = atol)
    end

    function Polyhedron(; verts = nothing, edges = nothing, facets = nothing)
        if isnothing(facets)
            facets = Vector{Int64}[]
        end

        if isnothing(edges)
            if facets == []
                edges = Vector{Int64}[]
            else
                edges = Vector{typeof(facets[1][1])}[]
                for f in facets
                    n = length(f)
                    append!(edges, [f[[mod1(i, n), mod1(i+1, n)]] for i in 1:n])
                end
                edges = collect.(unique(Set.(edges)))
            end
        end

        if isnothing(verts)
            error("verts is nothing. Case not implemented yet.")
        end

        return Polyhedron(verts, edges, facets)
    end
end


get_verts(poly::AbstractEmbOrCombPolyhedron) = deepcopy(poly.verts)

function set_verts!(poly::AbstractPolyhedron, verts::Matrix{<:Real})
    d = size(verts)[1]
    # @assert d == 3 "Only 3-dimensional polyhedra supported."
    poly.verts = verts
end

function set_verts!(poly::AbstractPolyhedron, verts::Vector{<:Vector{<:Real}})
    set_verts!(poly, hcat(verts...))
end

get_edges(poly::AbstractEmbOrCombPolyhedron) = deepcopy(poly.edges)

function set_edges!(poly::AbstractEmbOrCombPolyhedron, edges::Vector{<:Vector{<:Int}})
    @assert all(length(e) == 2 for e in edges) "Edges need to consist of vectors of length 2."
    # @assert sort(union(edges...)) == [1:max(union(edges...)...)...] "Vertex indices need to be 1, ..., $(length(unique(vcat(edges...))))."
    # TODO: Assert, dass die Kanten auch tatsächlich auf Rand von Facets liegen?
    poly.edges = edges
end

function get_facets(poly::AbstractEmbOrCombPolyhedron)
    return deepcopy(poly.facets)
end

function set_facets!(poly::AbstractEmbOrCombPolyhedron, facets::Vector{<:Vector{<:Int}}; atol::Real = 1e-8)
    # @assert sort(union(facets...)) == [1:max(union(facets...)...)...] "Vertex indices need to be 1, ..., $(length(unique(vcat(facets...))))."
    if any([affinedim(get_verts(poly)[:, f]; atol = atol) != 2 for f in facets])
        error("Facets have to span affine spaces of dimension $(2).")
    end
    for j = 1:length(facets)
        facet1 = facets[j]
        edges1 = Set.(union([[facet1[i], facet1[i+1]] for i in 1:length(facet1)-1], [[facet1[end], facet1[1]]]))
        for k = j+1:length(facets)
            facet2 = facets[k]
            edges2 = Set.(union([[facet2[i], facet2[i+1]] for i in 1:length(facet2)-1], [[facet2[end], facet2[1]]]))
            if setdiff(edges1, edges2) == [] && setdiff(edges2, edges1) == [] 
                error("One of facets $(facet1) and $(facet2) is contained in the other.")
            end
        end
    end
    poly.facets = facets
end

function ==(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron; atol = 1e-12)
    verts1 = get_verts(poly1)
    verts2 = get_verts(poly2)
    if size(verts1)[2] != size(verts2)[2]
        return false
    end

    edges1 = get_edges(poly1)
    edges2 = get_edges(poly2)
    if length(edges1) != length(edges2)
        return false
    end
    if length(Base.intersect(Set.(edges1), Set.(edges2))) < length(edges1)
        return false
    end

    facets1 = get_facets(poly1)
    facets2 = get_facets(poly2)
    if length(facets1) != length(facets2)
        return false
    end
    if length(Base.intersect(Set.(facets1), Set.(facets2))) < length(facets1)
        return false
    end
    
    return true
end

"""
    iscongruent(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron)

i.e. they have the same combinatorics and there exists a rigid map mapping the verts of poly1 to the verts of poly2.
"""
function iscongruent(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron)
    if Set.(get_edges(poly1)) != Set.(get_edges(poly2))
        return false
    end

    if !all([facet in get_facets(poly2) || reverse(facet) in get_facets(poly2) for facet in get_facets(poly1)])
        return false
    end

    try
        aff = rigidmap(get_verts(poly1), get_verts(poly2))
    catch AssertionError
        return false
    end

    return true
end

"""
    dimension(poly::AbstractPolyhedron)

Get the dimension of the unerlying space the polyhedron is embedded into.
"""
function dimension(poly::AbstractPolyhedron)
    return size(get_verts(poly))[1]
end


function display(poly::AbstractPolyhedron)
    print("""$(typeof(poly)) embedded into $(dimension(poly))-space with $(size(get_verts(poly))[2]) vertices, $(length(get_edges(poly))) edges and $(length(get_facets(poly))) facets.\n 
    Edges:  $(get_edges(poly))
    Facets: $(get_facets(poly)) \n""")
end

function display(poly::AbstractCombPolyhedron)
    print("""$(typeof(poly)) with $(length(get_verts(poly))) vertices, $(length(get_edges(poly))) edges and $(length(get_facets(poly))) facets.\n 
    Edges:  $(get_edges(poly))
    Facets: $(get_facets(poly)) \n""")
end


"""
    isadjacent(poly::AbstractEmbOrCombPolyhedron, facetoredge::AbstractVector{<:Integer}, facet::AbstractVector{<:Integer}; check::Bool = true)

Calculate whether the facet or edge facetoredge and the facet facet of the polyhedron poly are adjacent, i.e. share a common edge. 
If check is set to false, the function will not check whether facetoredge and facet are actually facets or edges of poly.
"""
function isadjacent(poly::AbstractEmbOrCombPolyhedron, facetoredge::AbstractVector{<:Integer}, facet::AbstractVector{<:Integer}; check::Bool = true)
    if check
        @assert Set(facetoredge) in Set.(get_edges(poly)) || Set(facetoredge) in Set.(get_facets(poly)) "facetoredge is not an edge or facet of poly."
        @assert Set(facet) in Set.(get_facets(poly)) "facet is not a facet of poly."
    end

    intersection = Base.intersect(facetoredge, facet)
    if length(intersection) < 2 || length(intersection) == length(facet)
        return false
    end

    ind_1 = indexin(intersection, facetoredge)
    ind_2 = indexin(intersection, facet)

    return (abs(ind_1[1] - ind_1[2]) == 1 || Set(facetoredge[ind_1]) == Set(facetoredge[[1,end]])) && (abs(ind_2[1] - ind_2[2]) == 1 || Set(facet[ind_2]) == Set(facet[[1,end]]))
end


"""
    adjfacets(poly::AbstractEmbOrCombPolyhedron, facetoredge::AbstractVector{<:Integer}; check::Bool = true)

Calculate the adjacent facets of the facet or edge facetoredge in the Polyhedron poly, i.e. the facets sharing at least one edge with facet.
If check is set to false, the function will not check whether facetoredge is actually a facet or edge of poly.
"""
function adjfacets(poly::AbstractEmbOrCombPolyhedron, facetoredge::AbstractVector{<:Integer}; check::Bool = true)
    if check
        @assert Set(facetoredge) in Set.(get_edges(poly)) || Set(facetoredge) in Set.(get_facets(poly)) "facetoredge is not an edge or facet of poly."
    end
    
    return filter(f -> isadjacent(poly, facetoredge, f, check = false), get_facets(poly))
end


"""
    isincident(v::Integer, facetoredge::AbstractVector{<:Integer})

Return true if the vertex with index v is incident to the facet or edge facetoredge, i.e. it is contained in facetoredge.
"""
function isincident(v::Integer, facetoredge::AbstractVector{<:Integer})
    return v in facetoredge
end


"""
    incfacets(poly::AbstractPolyhedron, vertexarray::AbstractVector{<:Integer})

Return the facets of the polyhedron poly, that are incident to all vertices in vertexarray, i.e. that contain all vertices in vertexarray.
"""
function incfacets(poly::AbstractEmbOrCombPolyhedron, vertexarray::AbstractVector{<:Integer})
    return filter(f -> all([isincident(v, f) for v in vertexarray]), get_facets(poly))
end


"""
    incfacets(poly::AbstractPolyhedron, v::Integer)

Return the facets of the polyhedron poly, that are incident to the vertex v, i.e. that contain v.
"""
function incfacets(poly::AbstractEmbOrCombPolyhedron, v::Integer)
    return filter(f -> isincident(v, f), get_facets(poly))
end


"""
    incedges(poly::AbstractPolyhedron, v::Integer)

Return the edges of the polyhedron poly, that are incident to the vertex v, i.e. that contain v.
"""
function incedges(poly::AbstractEmbOrCombPolyhedron, v::Integer)
    return filter(e -> isincident(v, e), get_edges(poly))
end

"""
    incedges(poly::AbstractPolyhedron, f::AbstractVector{<:Integer})

Return the edges of the polyhedron poly, that are incident to the facet f, i.e. that are a subset of f. If check is set to true, it is checked, whether f is a facet of poly.
"""
function incedges(poly::AbstractEmbOrCombPolyhedron, f::AbstractVector{<:Integer}; check::Bool = true)
    if check
        @assert Set(f) in Set.(get_facets(poly)) "f is not a facet of poly."
    end
    return filter(e -> isadjacent(poly, e, f, check = false), get_edges(poly))
end


function edge_direction(e::AbstractVector{<:Integer}, f::AbstractVector{<:Integer})
    @assert length(e) == 2 "e has to be a vector of length 2, but got $e."
    @assert issubset(e, f) "e ($e) is not an edge of f ($f)"
    ind = indexin(e, f)
    n = length(f)
    if ind == [n,1] || ind[1] + 1 == ind[2]
        return 1
    elseif ind == [1,n] || ind[1] - 1 == ind[2]
        return -1
    end

    error("e ($e) is not an edge of f ($f)")
end


"""
    boundary(poly::AbstractPolyhedron)

Compute the boundary of the polyhedron poly, i.e. the set of edges that are incident to exactly one facet.
"""
function boundary(poly::AbstractEmbOrCombPolyhedron)
    return filter(e -> length(adjfacets(poly, e)) == 1, get_edges(poly))
end

function removefacet!(poly::AbstractEmbOrCombPolyhedron, f::AbstractVector{<:Integer})
    @assert Set(f) in Set.(get_facets(poly)) "Facet f not in poly."
    fac = filter(facet -> Base.intersect(f, facet) == f, get_facets(poly))[1]
    relevantedges = incedges(poly, fac)
    setdiff!(poly.facets, [fac])
    setdiff!(poly.edges, Base.intersect(relevantedges, boundary(poly)))
end

function removefacet(poly::AbstractEmbOrCombPolyhedron, f::AbstractVector{<:Integer})
    p = deepcopy(poly)
    removefacet!(p, f)
    return p
end

function removeedge!(poly::AbstractEmbOrCombPolyhedron, e::AbstractVector{<:Integer})
    @assert Set(e) in Set.(get_edges(poly)) "Edge e not in poly."
    ed = filter(edge -> Base.intersect(e, edge) == e, get_edges(poly))[1]
    relevantfacets = incfacets(poly, ed)
    for f in relevantfacets
        removefacet!(poly, f)
    end
    setdiff!(poly.edges, [ed])
end

function removeedge(poly::AbstractEmbOrCombPolyhedron, e::AbstractVector{<:Integer})
    p = deepcopy(poly)
    removeedge!(p, e)
    return p
end

"""
    formpath!(vertexarray::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}})

Sort the vector vertexarray such that they lie on a common path along the edges defined in the vector edges.
"""
function formpath!(vertexarray::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}})
    endpoints = filter(v -> length(filter(e -> v in e, edges)) == 1, vertexarray)
    if length(endpoints) > 2
        @info "vertexarray: $(vertexarray)"
        @info "relevant edges: $(edges[[length(Base.intersect(e, vertexarray)) == 2 for e in edges]])"
        @info "endpoints: $(endpoints)"
    end
    @assert length(endpoints) <= 2 "Vertices don't lie on a common path"
    intersectionverts = filter(v -> length(filter(e -> v in e, edges)) > 2, vertexarray)
    if length(intersectionverts) > 0
        @info "vertexarray: $(vertexarray)"
        @info "relevant edges: $(edges[[length(Base.intersect(e, vertexarray)) == 2 for e in edges]])"
        @info "intersectionverts: $(intersectionverts)"
    end
    
    @assert length(intersectionverts) == 0 "No intersections allowed."

    tosort = deepcopy(vertexarray)
    if length(endpoints) == 0
        start = vertexarray[1]
    else
        start = endpoints[1]
    end
    path = [start]
    setdiff!(tosort, [start])
    vertexarray[1] = start

    for i = 2:length(vertexarray)
        next = tosort[map( j -> (Set([vertexarray[i-1], j]) in Set.(edges)), tosort )][1]
        vertexarray[i] = next
        setdiff!(tosort, [next])
    end
end


"""
    formpath(vertexarray::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}})

Sort the vector vertexarray such that they lie on a common path along the edges defined in the vector edges.
"""
function formpath(vertexarray::AbstractVector{<:Integer}, edges::AbstractVector{<:AbstractVector{<:Integer}})
    vertexarraycopy = deepcopy(vertexarray)
    formpath!(vertexarraycopy, edges)
    return vertexarraycopy
end

"""
    formpath!(vertexindices::Vector{<:Int}, poly::AbstractPolyhedron)

Sort the vector vertexindices such that the corresponding vertices of the Polyhedron poly form a vertex edges path.
"""
function formpath!(vertexindices::Vector{<:Int}, poly::AbstractEmbOrCombPolyhedron)
    edges = get_edges(poly)
    formpath!(vertexindices, edges)
end

"""
    formpath(vertexindices::AbstractVector{<:Integer}, poly::AbstractPolyhedron)

Sort the vector vertexindices such that the corresponding vertices of the Polyhedron poly form a vertex edges path.
"""
function formpath(vertexindices::AbstractVector{<:Integer}, poly::AbstractEmbOrCombPolyhedron)
    vertexindicescopy = deepcopy(vertexindices)
    formpath!(vertexindicescopy, poly)
    return vertexindicescopy
end


# # decide whether point is inside of polyhedron
# """
#     inpolyhedron(point::AbstractVector{<:Real}, poly::AbstractPolyhedron; atol::Real=1e-8)::Int


#     Randomized algorithm to check whether a point is contained in a polyhedron.
# """
# function inpolyhedron(point::AbstractVector{<:Real}, poly::AbstractPolyhedron; atol::Real=1e-8)::Int 
#     # check whether point lies on the boundary of poly
#     for facet in get_facets(poly)
#         polygon = get_verts(poly)[:, facet]
#         if  inpolygon3d(polygon, point, atol= atol) != 0
#             return -1
#         end
#     end

#     while true
#         v = normalize!(rand(Float64, 3))
#         # println(v)
#         r = Ray(point, v)
#         numIntersections = 0 # number of intersections of r and the facets of poly

#         for facet in get_facets(poly)
#             E = Plane(get_verts(poly)[:, facet])

#             try
#                 intersect(r, E)
#             catch error
#                 continue
#             end

#             p = intersect(r, E)

#             if inpolygon3d(get_verts(poly)[:, facet], p, atol = atol) == -1
#                 error("Ray intersects the boundary of a facet.")
#                 break
#             elseif inpolygon3d(get_verts(poly)[:, facet], p, atol = atol) == 1
#                 numIntersections = numIntersections + 1
#             end
#         end

#         # println(numIntersections)
#         if mod(numIntersections, 2) == 0
#             return 0
#         else
#             return 1
#         end
#     end
# end


"""
    orient_facets!(poly::AbstractPolyhedron; atol::Real = 1e-8)

Orient the facets of the polyhedron poly in line with the first facet. 
"""
function orient_facets!(poly::AbstractEmbOrCombPolyhedron; atol::Real = 1e-8)
    # https://stackoverflow.com/questions/48093451/calculating-outward-normal-of-a-non-convex-polyhedral#comment83177517_48093451
    adj = facet_adjacency(poly)

    oriented = [1]
    ext_facets = Dict{Int, Int}() # exterior facets that have to be oriented. key is facet index to be oriented, value is facet wrt which the key has to be oriented.
    for i in findall(adj[:, 1])
        ext_facets[i] = 1
    end

    while length(oriented) < length(get_facets(poly))
        f_ind = collect(keys(ext_facets))[1]
        f = get_facets(poly)[f_ind]
        g_ind = ext_facets[f_ind]
        g = get_facets(poly)[ext_facets[f_ind]]
        
        inter = Base.intersect(f,g)
        # indin_f = indexin(inter, f)
        # indin_g = indexin(inter, g)
        
        # delete f_ind from ext_facets and add all adjacent facets of f to ext_facets
        delete!(ext_facets, f_ind)
        for i in Base.intersect(findall(adj[:, f_ind]), setdiff(collect(1:length(get_facets(poly))), oriented))
            ext_facets[i] = f_ind
        end
        
        # add f_ind to oriented
        oriented = [oriented; f_ind]
        
        # orient f with regard to g
        if edge_direction(inter, f) == edge_direction(inter, g)
            poly.facets[f_ind] = reverse(poly.facets[f_ind])
        end
    end

    return poly
end


"""
    orient_facets(poly::AbstractPolyhedron)

Orient the facets of the polyhedron poly.
"""
function orient_facets(poly::AbstractEmbOrCombPolyhedron; atol::Real = 1e-8)
    polycopy = deepcopy(poly)
    orient_facets!(polycopy; atol = atol)
    return polycopy
end

"""
    signed_vol(poly::AbstractPolyhedron; atol::Real = 1e-8, is_oriented = false)

Calculate the signed volume of the polyhedron according to the orientation of the facets.
If is_oriented = false, calculate an orientation of the polyhedron first.
"""
function vol_signed(poly::AbstractPolyhedron; atol::Real = 1e-8, is_oriented = false)
    if dimension(poly) != 2 && dimension(poly) != 3
        error("Volume only defined in R² and R³.")
    end

    if !is_oriented
        poly_orient = orient_facets(poly)
    else 
        poly_orient = deepcopy(poly)
    end

    # now all facets are oriented. Volume can be calculated from the signed volumes of parallelopiped
    vol = 0
    if dimension(poly_orient) == 2
        # https://en.wikipedia.org/wiki/Polygon#Simple_polygons
        for f in get_facets(poly_orient)
            n = length(f)
            for k in 1:length(f)
                # signed volume of the triangle spanned by [0,0,0], f[1], f[k], f[k+1] (=1/2*volume of parallelogram). The sum of all those for every facet is the signed vol of poly.
                vol += 1/2 * det(get_verts(poly_orient)[:, [f[mod1(k,n)], f[mod1(k+1, n)]]])
            end
        end
    elseif dimension(poly_orient) == 3
        for f in get_facets(poly_orient)
            for k in 2:length(f)-1
                # signed volume of the tetrahedron spanned by [0,0,0], f[1], f[k], f[k+1] (=1/6*volume of parallelopiped). The sum of all those for every facet is the signed vol of poly.
                vol += 1/6 * det(get_verts(poly_orient)[:, [f[1], f[k], f[k+1]]])
            end
        end
    end

    return vol
end

"""
    orient_facets_outwardn!(poly::AbstractPolyhedron; atol::Real = 1e-8)

Orient the facets so that the facets are oriented counterclockwise wrt the outward normals.
If is_oriented == true, the poly is expected to be oriented.
"""
function orient_facets_ccw!(poly::AbstractPolyhedron; atol::Real = 1e-8, is_oriented::Bool = false)
    if !is_oriented
        poly_orient = orient_facets(poly)
    else 
        poly_orient = poly
    end

    vol = vol_signed(poly_orient)
    # if sign is negative, the facets are ordered clockwise wrt the outward normal
    if vol >= 0
        set_facets!(poly, get_facets(poly_orient); atol = atol)
    else
        set_facets!(poly, reverse.(get_facets(poly_orient)); atol = atol)
    end
end

function orient_facets_ccw(poly::AbstractPolyhedron; atol::Real = 1e-8, is_oriented::Bool = false)
    poly_copy = deepcopy(poly)
    orient_facets_ccw!(poly_copy; atol = atol, is_oriented = is_oriented)
    return poly_copy
end

"""
    vol(poly::AbstractPolyhedron)

Calculate the volume of the polyhedron poly.
If is_oriented == true, poly is expected to be oriented.
"""
function vol(poly::AbstractPolyhedron; atol::Real = 1e-8, is_oriented::Bool = false)
    return abs(vol_signed(poly, atol = atol, is_oriented = is_oriented))
end

"""
    facet_adjacency(poly::AbstractPolyhedron)

Adjacency matrix for the facets of the polyhedron poly. Two facets are adjacent if they share an edge.
"""
function facet_adjacency(poly::AbstractEmbOrCombPolyhedron)
    facets = get_facets(poly)
    adj = zeros(Bool, length(facets), length(facets))
    for i in 1:length(facets)
        for j in i+1:length(facets)
            if isadjacent(poly, facets[i], facets[j], check = false)
                adj[i,j] = true
                adj[j,i] = true
            end
        end
    end
    return adj
end