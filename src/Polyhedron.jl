using GenericLinearAlgebra
using PolygonOps
import Base.Multimedia.display
import Base.==

include("Framework.jl")
include("affine_geometry.jl")
include("polygonal_geometry.jl")

abstract type AbstractPolyhedron{S<:Real,T<:Integer} <: AbstractEmbeddedGraph{S, T} end

mutable struct Polyhedron{S<:Real, T<:Integer} <:AbstractPolyhedron{S, T}
    verts::Matrix{S} # matrix of coordinates of the vertices. Columns correspond to vertices.
    edges::Vector{Vector{T}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{Vector{T}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
    # TODO: Neuen constructor mit optionalen Argumenten (wenn nur coordinates gegeben werden, ist Ergebnis die konvexe Hülle der Punkte + check der Dimension
    # wenn nur Facets gegeben sind, werden Edges automatisch gesetzt und es wird gecheckt, dass Vertizes auf einer Facet koplanar aber nicht kollinear sind)
    function Polyhedron(verts::AbstractMatrix{<:Real}, edges::Vector{<:Vector{<:Integer}}, facets::Vector{<:Vector{<:Integer}}; atol::Real = 1e-8)
        framework = Framework(verts, edges)
        if any([affinedim(verts[:,f]; atol = atol) != 2 for f in facets])
            error("Facets have to span a space of affine dimension 2.")
        end
        # TODO: Facets sind zyklische Graphen -> In Framework.jl für Graphen implementieren: iscyclic.
        S = typeof(verts[1,1])
        T = typeof(edges[1][1])
        return orient_facets(new{S,T}(verts, edges, facets); atol = atol)
    end

    function Polyhedron(verts::Vector{<:Vector{<:Real}}, edges::Vector{<:Vector{<:Integer}}, facets::Vector{<:Vector{<:Integer}}; atol::Real = 1e-8)
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
                    append!(edges, [[f[mod1(i, n)], f[mod1(i+1, n)]] for i in 1:n])
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


function get_verts(poly::AbstractPolyhedron)
    return deepcopy(poly.verts)
end

function set_verts!(poly::AbstractPolyhedron, verts::Matrix{<:Real})
    d = size(verts)[1]
    # @assert d == 3 "Only 3-dimensional polyhedra supported."
    poly.verts = verts
end

function set_verts!(poly::AbstractPolyhedron, verts::Vector{<:Vector{<:Real}})
    set_verts!(poly, hcat(verts...))
end


function set_edges!(poly::AbstractPolyhedron, edges::Vector{<:Vector{<:Int}})
    @assert all(length(e) == 2 for e in edges) "Edges need to consist of vectors of length 2."
    # @assert sort(union(edges...)) == [1:max(union(edges...)...)...] "Vertex indices need to be 1, ..., $(length(unique(vcat(edges...))))."
    # TODO: Assert, dass die Kanten auch tatsächlich auf Rand von Facets liegen?
    poly.edges = edges
end


function get_facets(poly::AbstractPolyhedron)
    return deepcopy(poly.facets)
end

function set_facets!(poly::AbstractPolyhedron, facets::Vector{<:Vector{<:Int}}; atol::Real = 1e-8)
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

function remove_edge!(poly::AbstractPolyhedron, e::Vector{<:Int}; atol = 1e-8)
    if e in get_edges(poly)
        edge = e
    elseif reverse(e) in get_edges(poly)
        edge = reverse(e)
    else
        error("Edge not in polyhedron.")
    end
    neighbors = adjfacets(poly, edge)

    if length(neighbors) == 1
        # edge can be removed by removing the whole adjacent facet
        border_edges = filter(e -> length(adjfacets(poly, e)) == 1, incedges(poly, neighbors[1]))
        setdiff!(poly.facets, neighbors)
        setdiff!(poly.edges, border_edges)
        # border edges describe a vertex edge path between the two neighbors. All inner vertices of the path need to be removed. Start and end point remain.
        endpoints = filter(v -> count(x -> x == v, vcat(border_edges...)) == 1, vcat(border_edges...))
        remove_verts = unique(setdiff(vcat(border_edges...), endpoints))
    elseif length(neighbors) == 2
        # edge can only be removed, if neighboring facets are coplanar
        if affinedim(get_verts(poly)[:, unique(vcat(neighbors...))], atol = atol) != 2
            error("Edge can only be removed if neighboring facets are coplanar.")
        end

        # if one edge between neighbors is removed, all edges between neighbors are removed
        border_edges = Base.intersect(incedges(poly, neighbors[1]), incedges(poly, neighbors[2]))

        # border edges describe a vertex edge path between the two neighbors. All inner vertices of the path need to be removed. Start and end point remain.
        endpoints = filter(v -> count(x -> x == v, vcat(border_edges...)) == 1, vcat(border_edges...))
        start = endpoints[1]
        finish = endpoints[2]
        remove_verts = unique(setdiff(vcat(border_edges...), endpoints))

        # remove the vertices to be removed from neighbors. Combine them together so that the new facet is the union of the two old facets withouth the removed vertices.
        setdiff!(poly.facets, neighbors) # remove neighbors from facets to add them combined later
        neighbors = [setdiff(neighbors[1], remove_verts), setdiff(neighbors[2], remove_verts)]
        # rearrange neighbors so that in neighbors[1] the start vertex is the first and the finish vertex is the last and in neighbors[2] it's the other way around.
        start_index1 = findfirst(x -> x == start, neighbors[1])
        finish_index2 = findfirst(x -> x == finish, neighbors[2])

        neighbors[1] = neighbors[1][[mod1(i+start_index1-1, length(neighbors[1])) for i in eachindex(neighbors[1])]] # now first entry of neighbors[1] is start
        if neighbors[1][end] != finish # if last entry of neighbors[1] is not finish, then finish is the second entry of neighbors[1]. Reverse neighbors[1] and shift by one to the right so that first entry is start and last one is finish.
            neighbors[1] = reverse(neighbors[1])[[mod1(i-1, length(neighbors[1])) for i in eachindex(neighbors[1])]]
        end

        neighbors[2] = neighbors[2][[mod1(i+finish_index2-1, length(neighbors[2])) for i in eachindex(neighbors[2])]] # now first entry of neighbors[2] is finish
        if neighbors[2][end] != start # if last entry of neighbors[2] is not start, then start is the second entry of neighbors[2]. Reverse neighbors[2] and shift by one to the right so that first entry is finish and last one is start.
            neighbors[2] = reverse(neighbors[2])[[mod1(i-1, length(neighbors[2])) for i in eachindex(neighbors[2])]]
        end

        # now the vertex orders of neighbors[1] and neighbors[2] are consistent with the vertex order of the new facet.
        newfacet = unique(vcat(neighbors[1], neighbors[2]))
        
        # add new facet to poly
        push!(poly.facets, newfacet)
    else
        error("Edge has more than 2 neighboring facets.")
    end

    # remove the border edges from poly
    setdiff!(poly.edges, border_edges)

    # remove remove_verts from poly by deleting the columns from verts and shifting the corresponding entries in the edges and facets so that the vertex labels are 1,...,n
    vertexmap(v) = v - count(x -> x < v, remove_verts)
    poly.edges = [vertexmap.(e) for e in poly.edges]
    poly.facets = [vertexmap.(f) for f in poly.facets]
    poly.verts = poly.verts[:, setdiff(1:size(poly.verts)[2], remove_verts)]

    return poly
end

function remove_edge(poly::AbstractPolyhedron, e::Vector{<:Int}; atol = 1e-8)
    p = deepcopy(poly)
    remove_edge!(p, e; atol = atol)
    return p
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


"""
    isadjacent(poly::AbstractPolyhedron, facetoredge::AbstractVector{<:Integer}, facet::AbstractVector{<:Integer})

Calculate whether the facet or edge facetoredge and the facet facet of the polyhedron poly are adjacent, i.e. share a common edge.
"""
function isadjacent(poly::AbstractPolyhedron, facetoredge::AbstractVector{<:Integer}, facet::AbstractVector{<:Integer})
    @assert Set(facetoredge) in Set.(get_facets(poly)) || Set(facetoredge) in Set.(get_edges(poly)) "facetoredge has to be a facet or an edge of poly."
    @assert Set(facet) in Set.(get_facets(poly)) "facet has to be a facet of poly."

    intersection = Base.intersect(facetoredge, facet)

    return any(map(edge -> Base.intersect(edge, intersection) == edge, get_edges(poly)))
end


"""
    adjfacets(poly::AbstractPolyhedron, facetoredge::AbstractVector{<:Integer})

Calculate the adjacent facets of the facet or edge facetoredge in the Polyhedron poly, i.e. the facets sharing at least one edge with facet.
"""
function adjfacets(poly::AbstractPolyhedron, facetoredge::AbstractVector{<:Integer})
    sol = filter(facet2 -> isadjacent(poly, facetoredge, facet2) && Set(facet2) != Set(facetoredge), get_facets(poly))
    setdiff!(sol, [facetoredge])
    return sol
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
function incfacets(poly::AbstractPolyhedron, vertexarray::AbstractVector{<:Integer})
    return filter(f -> all([isincident(v, f) for v in vertexarray]), get_facets(poly))
end


"""
    incfacets(poly::AbstractPolyhedron, v::Integer)

Return the facets of the polyhedron poly, that are incident to the vertex v, i.e. that contain v.
"""
function incfacets(poly::AbstractPolyhedron, v::Integer)
    return filter(f -> isincident(v, f), get_facets(poly))
end


"""
    incedges(poly::AbstractPolyhedron, v::Integer)

Return the edges of the polyhedron poly, that are incident to the vertex v, i.e. that contain v.
"""
function incedges(poly::AbstractPolyhedron, v::Integer)
    return filter(e -> isincident(v, e), get_edges(poly))
end

"""
    incedges(poly::AbstractPolyhedron, f::AbstractVector{<:Integer})

Return the edges of the polyhedron poly, that are incident to the facet f, i.e. that are a subset of f.
"""
function incedges(poly::AbstractPolyhedron, f::AbstractVector{<:Integer})
    return filter(e -> issubset(e, f), get_edges(poly))
end


"""
    boundary(poly::AbstractPolyhedron)

Compute the boundary of the polyhedron poly, i.e. the set of facets that are incident to exactly one edge.
"""
function boundary(poly::AbstractPolyhedron)
    return filter(f -> length(incedges(poly, f)) == 1, get_facets(poly))
end

function removefacet!(poly::AbstractPolyhedron, f::AbstractVector{<:Integer})
    @assert Set(f) in Set.(get_facets(poly)) "Facet f not in poly."
    fac = filter(facet -> Base.intersect(f, facet) == f, get_facets(poly))[1]
    relevantedges = incedges(poly, fac)
    setdiff!(poly.facets, [fac])
    setdiff!(poly.edges, Base.intersect(relevantedges, boundary(poly)))
end

function removefacet(poly::AbstractPolyhedron, f::AbstractVector{<:Integer})
    p = deepcopy(poly)
    removefacet!(p, f)
    return p
end

function removeedge!(poly::AbstractPolyhedron, e::AbstractVector{<:Integer})
    @assert Set(e) in Set.(get_edges(poly)) "Edge e not in poly."
    ed = filter(edge -> Base.intersect(e, edge) == e, get_edges(poly))[1]
    relevantfacets = incfacets(poly, ed)
    for f in relevantfacets
        removefacet!(poly, f)
    end
    setdiff!(poly.edges, [ed])
end

function removeedge(poly::AbstractPolyhedron, e::AbstractVector{<:Integer})
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
function formpath!(vertexindices::Vector{<:Int}, poly::AbstractPolyhedron)
    edges = get_edges(poly)
    formpath!(vertexindices, edges)
end

"""
    formpath(vertexindices::AbstractVector{<:Integer}, poly::AbstractPolyhedron)

Sort the vector vertexindices such that the corresponding vertices of the Polyhedron poly form a vertex edges path.
"""
function formpath(vertexindices::AbstractVector{<:Integer}, poly::AbstractPolyhedron)
    vertexindicescopy = deepcopy(vertexindices)
    formpath!(vertexindicescopy, poly)
    return vertexindicescopy
end


# decide whether point is inside of polyhedron
"""
    inpolyhedron(point::AbstractVector{<:Real}, poly::AbstractPolyhedron; atol::Real=1e-8)::Int


    Randomized algorithm to check whether a point is contained in a polyhedron.
"""
function inpolyhedron(point::AbstractVector{<:Real}, poly::AbstractPolyhedron; atol::Real=1e-8)::Int 
    # check whether point lies on the boundary of poly
    for facet in get_facets(poly)
        polygon = get_verts(poly)[:, facet]
        if  inpolygon3d(polygon, point, atol= atol) != 0
            return -1
        end
    end

    while true
        v = normalize!(rand(Float64, 3))
        # println(v)
        r = Ray(point, v)
        numIntersections = 0 # number of intersections of r and the facets of poly

        for facet in get_facets(poly)
            E = Plane(get_verts(poly)[:, facet])

            try
                intersect(r, E)
            catch error
                continue
            end

            p = intersect(r, E)

            if inpolygon3d(get_verts(poly)[:, facet], p, atol = atol) == -1
                error("Ray intersects the boundary of a facet.")
                break
            elseif inpolygon3d(get_verts(poly)[:, facet], p, atol = atol) == 1
                numIntersections = numIntersections + 1
            end
        end

        # println(numIntersections)
        if mod(numIntersections, 2) == 0
            return 0
        else
            return 1
        end
    end
end


"""
    orient_facets!(poly::AbstractPolyhedron)

Orient the facets of the polyhedron poly.
"""
function orient_facets!(poly::AbstractPolyhedron; atol::Real = 1e-8)
    # https://stackoverflow.com/questions/48093451/calculating-outward-normal-of-a-non-convex-polyhedral#comment83177517_48093451
    newfacets = Vector{Int}[]

    f = get_facets(poly)[1]
    adj = adjfacets(poly, f)
    push!(newfacets, f)

    # exterior facets that have been oriented
    extfacets = Vector{Int}[]

    function direction(facet, edge)
        # function that determines if edge is oriented forwards (1) or backwards (-1) in facet. 
        inds = indexin(edge, facet)
        if mod1(inds[2] - inds[1], length(facet)) == 1 return 1 end
        return -1
    end

    while length(newfacets) < length(get_facets(poly))
        setdiff!(extfacets, [f])

        for g in adj
            if g in newfacets || reverse(g) in newfacets 
                continue
            end

            inter = Base.intersect(f,g)
            e = [inter[1], inter[2]]

            if direction(f, [inter[1], inter[2]]) == direction(g, [inter[1], inter[2]])
                reverse!(g)
            end

            # g is now oriented and a newfacet and an exteriorfacet in this step
            push!(newfacets, g)
            push!(extfacets, g)
        end

        for g in extfacets
            adj_new = adjfacets(poly, g)
            if length(setdiff(Set.(adj_new), Set.(newfacets))) > 0
                f = g
                adj = adj_new
                break
            end
            # there are no adjacent facets of g that are not oriented already -> g is not exterior.
            setdiff!(extfacets, [g])
        end
    end

    set_facets!(poly, newfacets; atol = atol)
end

"""
    orient_facets(poly::AbstractPolyhedron)

Orient the facets of the polyhedron poly.
"""
function orient_facets(poly::AbstractPolyhedron; atol::Real = 1e-8)
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
        poly_orient = deeocopy(poly)
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
    return abs(vol_signed(poly))
end