using GenericLinearAlgebra
using PolygonOps
import Base.Multimedia.display
import Base.==

include("affine_geometry.jl")
include("polygonal_geometry.jl")

abstract type AbstractPolyhedron end
mutable struct Polyhedron <:AbstractPolyhedron
    verts::Vector{<:Vector{<:Real}} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{<:Vector{<:Integer}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{<:Vector{<:Integer}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
    # TODO: Neuen constructor mit optionalen Argumenten (wenn nur coordinates gegeben werden, ist Ergebnis die konvexe Hülle der Punkte + check der Dimension
    # wenn nur Facets gegeben sind, werden Edges automatisch gesetzt und es wird gecheckt, dass Vertizes auf einer Facet koplanar aber nicht kollinear sind)
end

function get_verts(poly::Polyhedron)
    return deepcopy(poly.verts)
end

function set_verts!(poly::Polyhedron, verts::Vector{<:Vector{<:Real}})
    d = length(verts[1])
    @assert d == 3 "Only 3-dimensional polyhedra supported."
    @assert all([length(v) == d for v in verts]) "Dimension mismatch in vertices."
    poly.verts = verts
end

function get_edges(poly::Polyhedron)
    return deepcopy(poly.edges)
end

function set_edges!(poly::Polyhedron, edges::Vector{<:Vector{<:Int}})
    @assert all(length(e) == 2 for e in edges) "Edges need to consist of vectors of length 2."
    # @assert sort(union(edges...)) == [1:max(union(edges...)...)...] "Vertex indices need to be 1, ..., $(length(unique(vcat(edges...))))."
    # TODO: Assert, dass die Kanten auch tatsächlich auf Rand von Facets liegen?
    poly.edges = edges
end


function get_facets(poly::Polyhedron)
    return deepcopy(poly.facets)
end

function set_facets!(poly::Polyhedron, facets::Vector{<:Vector{<:Int}})
    # @assert sort(union(facets...)) == [1:max(union(facets...)...)...] "Vertex indices need to be 1, ..., $(length(unique(vcat(facets...))))."
    @assert all([affinedim(get_verts(poly)[f]) == 2 for f in facets]) "Facets have to span affine spaces of dimension 2."
    for j = 1:length(facets)
        facet1 = facets[j]
        edges1 = Set.(union([[facet1[i], facet1[i+1]] for i in 1:length(facet1)-1], [[facet1[end], facet1[1]]]))
        for k = j+1:length(facets)
            facet2 = facets[k]
            edges2 = Set.(union([[facet2[i], facet2[i+1]] for i in 1:length(facet2)-1], [[facet2[end], facet2[1]]]))
            @assert setdiff(edges1, edges2) != [] && setdiff(edges2, edges1) != [] "One of facets $(facet1) and $(facet2) is contained in the other."
        end
    end
    poly.facets = facets
end

function ==(poly1::Polyhedron, poly2::Polyhedron; atol = 1e-12)
    verts1 = get_verts(poly1)
    verts2 = get_verts(poly2)
    if length(verts1) != length(verts2)
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
    iscongruent(poly1::Polyhedron, poly2::Polyhedron)

Determine whether two polyhedra are congruent, 
i.e. they have the same combinatorics and there exists a rigid map mapping the verts of poly1 to the verts of poly2.
"""
function iscongruent(poly1::Polyhedron, poly2::Polyhedron)
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
    dimension(poly::Polyhedron)

Get the dimension of the unerlying space the polyhedron is embedded into.
"""
function dimension(poly::Polyhedron)
    return length(get_verts(poly)[1])
end


function display(poly::Polyhedron)
    print("""Polyhedron embedded into $(dimension(poly))-space with $(length(get_verts(poly))) vertices, $(length(get_edges(poly))) edges and $(length(get_facets(poly))) facets.\n 
    Edges:  $(get_edges(poly))
    Facets: $(get_facets(poly)) \n""")
end


"""
    isadjacent(poly::Polyhedron, facetoredge::Vector{<:Int}, facet::Vector{<:Int})


Calculate whether the facet or edge facetoredge and the facet facet of the polyhedron poly are adjacent, i.e. share a common edge.
"""
function isadjacent(poly::Polyhedron, facetoredge::Vector{<:Int}, facet::Vector{<:Int})
    @assert Set(facetoredge) in Set.(get_facets(poly)) || Set(facetoredge) in Set.(get_edges(poly)) "facetoredge has to be a facet or an edge of poly."
    @assert Set(facet) in Set.(get_facets(poly)) "facet2 has to be a facet of poly."

    intersection = Base.intersect(facetoredge, facet)

    return any(map(edge -> Base.intersect(edge, intersection) == edge, get_edges(poly)))
end


"""
    adjfacets(poly::Polyhedron, facetoredge::Vector{<:Int})

Calculate the adjacent facets of the facet or edge facetoredge in the Polyhedron poly, i.e. the facets sharing at least one edge with facet.
"""
function adjfacets(poly::Polyhedron, facetoredge::Vector{<:Int})
    sol = filter(facet2 -> isadjacent(poly, facetoredge, facet2), get_facets(poly))
    setdiff!(sol, [facetoredge])
    return sol
end


"""
    isincident(v::Int, facetoredge::Vector{<:Int})

Return true if the vertex with index v is incident to the facet or edge facetoredge, i.e. it is contained in facetoredge.
"""
function isincident(v::Int, facetoredge::Vector{<:Int})
    return v in facetoredge
end


"""
    incfacets(poly::Polyhedron, vertexarray::Vector{<:Int})

Return the facets of the polyhedron poly, that are incident to all vertices in vertexarray, i.e. that contain all vertices in vertexarray.
"""
function incfacets(poly::Polyhedron, vertexarray::Vector{<:Int})
    return filter(f -> all([isincident(v, f) for v in vertexarray]), get_facets(poly))
end


"""
    incfacets(poly::Polyhedron, v::Int)

Return the facets of the polyhedron poly, that are incident to the vertex v, i.e. that contain v.
"""
function incfacets(poly::Polyhedron, v::Int)
    return filter(f -> isincident(v, f), get_facets(poly))
end


"""
    formpath!(vertexarray::Vector{<:Int}, edges::Vector{<:Vector{<:Int}})

Sort the vector vertexarray such that they lie on a common path along the edges defined in the vector edges.
"""
function formpath!(vertexarray::Vector{<:Int}, edges::Vector{<:Vector{<:Int}})
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
    formpath(vertexarray::Vector{<:Int}, edges::Vector{<:Vector{<:Int}})

Sort the vector vertexarray such that they lie on a common path along the edges defined in the vector edges.
"""
function formpath(vertexarray::Vector{<:Int}, edges::Vector{<:Vector{<:Int}})
    vertexarraycopy = deepcopy(vertexarray)
    formpath!(vertexarraycopy, edges)
    return vertexarraycopy
end

"""
    formpath!(vertexindices::Vector{<:Int}, poly::Polyhedron)

Sort the vector vertexindices such that the corresponding vertices of the Polyhedron poly form a vertex edges path.
"""
function formpath!(vertexindices::Vector{<:Int}, poly::Polyhedron)
    edges = get_edges(poly)
    formpath!(vertexindices, edges)
end

"""
    formpath(vertexindices::Vector{<:Int}, poly::Polyhedron)

Sort the vector vertexindices such that the corresponding vertices of the Polyhedron poly form a vertex edges path.
"""
function formpath(vertexindices::Vector{<:Int}, poly::Polyhedron)
    vertexindicescopy = deepcopy(vertexindices)
    formpath!(vertexindicescopy, poly)
    return vertexindicescopy
end


# decide whether point is inside of polyhedron
"""
    Randomized algorithm to check whether a point is contained in a polyhedron.
"""
function inpolyhedron(point::Vector{<:Real}, poly::Polyhedron;  atol::Real=1e-5)::Int 
    # check whether point lies on the boundary of poly
    for facet in get_facets(poly)
        polygon = push!(map(v -> poly.verts[v], facet), poly.verts[facet[1]])
        if  inpolygon3d(point, polygon,  atol= atol) != 0
            return -1
        end
    end

    while true
        r = Ray(point, normalize!(randn(Float64, 3)))
        numIntersections = 0 # number of intersections of r and the facets of poly

        for facet in get_facets(poly)
            polygon = push!(map(v -> poly.verts[v], facet), poly.verts[facet[1]])
            E = Plane(polygon)

            try
                intersect(r, E)
            catch error
                continue
            end

            p = intersect(r, E)

            if inpolygon3d(p, polygon,  atol =  atol) == -1
                break
            elseif inpolygon3d(p, polygon,  atol =  atol) == 1
                numIntersections = numIntersections + 1
            end
        end

        if mod(numIntersections, 2) == 0
            return 0
        else
            return 1
        end
    end
end



