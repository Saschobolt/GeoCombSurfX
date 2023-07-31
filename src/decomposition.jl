include("Polyhedron.jl")
include("combinatorics.jl")


# triangulate surface of Polyhedron
"""
returns a polyhedron containing the vertices and edges of poly such that every facet is triangular.
"""
function triangulate!(poly::Polyhedron; atol = 1e-5)
    newVerts = deepcopy(get_verts(poly))
    newEdges = deepcopy(get_edges(poly))
    newFacets = Vector{Int}[]

    # triangulate every facet of poly by the earcut algorithm
    for facet in get_facets(poly)
        subfacet = deepcopy(facet)
        polygon = push!(map(i -> get_verts(poly)[i], subfacet), get_verts(poly)[subfacet[1]])

        while length(subfacet) > 3
            # try to cut triangle between 2nd and last vertex
            if inpolygon3d((get_verts(poly)[subfacet[end]] + get_verts(poly)[subfacet[2]]) / 2, polygon, tol = atol) == 1
                push!(newEdges, [subfacet[end], subfacet[2]])
                push!(newFacets, [subfacet[end], subfacet[1], subfacet[2]])
                subfacet = subfacet[2:end]
                continue
            # try to cut triangle between 1st and second to last vertex
            elseif inpolygon3d((get_verts(poly)[subfacet[end - 1]] + get_verts(poly)[subfacet[1]]) / 2, polygon, tol = atol) == 1
                push!(newEdges, [subfacet[end - 1], subfacet[1]])
                push!(newFacets, [subfacet[end - 1], subfacet[length(subfacet)], subfacet[1]])
                subfacet = subfacet[1:end - 1]
                continue
            else
                # every other pair of vertices that are 1 vertex apart
                for ind in 2:length(subfacet) - 1
                    if inpolygon3d((get_verts(poly)[subfacet[ind - 1]] + get_verts(poly)[subfacet[ind + 1]]) / 2, polygon, tol = atol) == 1
                        push!(newEdges, [subfacet[ind - 1], subfacet[ind + 1]])
                        push!(newFacets, [subfacet[ind - 1], subfacet[ind], subfacet[ind + 1]])
                        subfacet = subfacet[1:end .!= ind]
                        break
                    end
                end
            end
        end

        push!(newFacets, subfacet)
    end

    set_verts!(poly, newVerts)
    set_edges!(poly, newEdges)
    set_facets!(poly, newFacets)
end


function triangulate(poly::Polyhedron; atol = 1e-5)::Polyhedron
    polycopy = deepcopy(poly)
    triangulate!(polycopy, atol = atol)

    return polycopy
end


"""
    flattenfacets!(poly::Polyhedron; atol = 1e-5)::Polyhedron

Remove flat edges of the Polyhedron poly.
"""
function flattenfacets!(poly::Polyhedron; atol = 1e-5)
    edgeschecked = []
    
    # @info "edgestocheck: $(edgestocheck)"
    # determine flat edges by calculating the affine dimension of the union of two facets.
    for edge in Set.(get_edges(poly))
        if edge in edgeschecked
            continue
        end

        @assert length(adjfacets(poly, collect(edge))) <= 2 "Edge $(collect(edge)) is degenerate as it is edge of the facets $(adjfacets(poly, collect(edge)))"

        if length(adjfacets(poly, collect(edge))) < 2
            continue
        end

        facet1 = adjfacets(poly, collect(edge))[1]
        facet2 = adjfacets(poly, collect(edge))[2]
        # possible edges in the intersection of the facets. All of those are checked at once since the affine dimension of a facet is 2. 
        # So one of them is flat <=> all of them are flat
        intersect_verts = Base.intersect(facet1, facet2)
        intersect_edges = Base.intersect(Set.(EdgesOfFace(poly, facet1)), Set.(EdgesOfFace(poly, facet2)))
        
        append!(edgeschecked, intersect_edges)
        unique!(edgeschecked)

        d = affinedim(get_verts(poly)[union(facet1, facet2)])

        # if facets span a space of dimension 2 the edges in the intersection can be removed
        # the facets need to be merged
        if d == 2
            # mark edges in the intersection as checked
            append!(edgeschecked, intersect_edges)


            # remove flat edges from poly
            oldedges = get_edges(poly)
            newedges = collect.(setdiff(Set.(oldedges), intersect_edges))
            set_edges!(poly, newedges)

            # merge facets 
            oldfacets = get_facets(poly)
            newfacets = setdiff(oldfacets, [facet1, facet2])
            push!(newfacets, formpath(union(facet1, facet2), poly))
            set_facets!(poly, newfacets)
        end

        append!(edgeschecked, intersect_edges)
    end

    # if a vertex lies inside a polygon forming a facet, it can be ignored. 
    # Thus the entries of sol_edges have to be transformed in order to be labeled 1,...,n
    function indexmap(i::Integer)
        return indexin(i, sort(union(get_edges(poly)...)))[1]
    end

    # if a vertex lies inside a polygon forming a facet, it can be ignored
    set_verts!(poly, get_verts(poly)[sort(union(get_edges(poly)...))])    

    set_edges!(poly, [indexmap.(e) for e in get_edges(poly)])
    set_facets!(poly, [indexmap.(f) for f in get_facets(poly)])
end


function flattenfacets(poly::Polyhedron, atol = 1e-5)
    polycopy = deepcopy(poly)
    flattenfacets!(polycopy)
    return polycopy
end

"""
e::Vector{Int}
poly::Polyhedron
tol::Float64
Checks whether the edge e is turnable edge of poly. I.e. poly is 
    - a spherical simplicial surface and 
    - the line connecting the wingtips of the butterfly with inner edge e is contained in poly.
Floats with abs value <tol are considered zeros
"""
function isturnable(e::Vector{<:Int}, polyhedron::Polyhedron; tol::Real=1e-5)::Bool
    @assert all(length.(get_facets(polyhedron)).==3) "poly may only have triangle facets"
    @assert in(e, get_edges(polyhedron)) || in(reverse(e), get_edges(poly)) "e has to be an edge of poly"
    
    # adjacent facets to e
    butterfly = get_facets(polyhedron)[map(f -> length(Base.intersect(f, e)) == 2, get_facets(polyhedron))]
    # turned edge
    butterflytips = setdiff(union(butterfly...), e)
    mid = (get_verts(polyhedron)[butterflytips[1]] + get_verts(polyhedron)[butterflytips[2]]) / 2

    # if midpoint of turned edge is inside the polyhedron the edge is turnable
    if inpolyhedron(mid, flattenfacets(polyhedron), tol = tol) == 0
        return 0
    end

    return 1
end


"""
poly::Polyhedron
tol::Float64
returns whether poly is convex
Floats with abs value <tol are considered zeros
"""
function isconvex(poly::Polyhedron; tol::Real=1e-5)::Bool
    polyhedron = deepcopy(poly)

    # if poly is a non degenerate tetrahedron, it is convex
    if length(get_verts(polyhedron)) == 4 && length(get_edges(polyhedron)) == 6 && length(unique(get_verts(polyhedron))) == length(get_verts(polyhedron))
        return 1
    end

    # otherwise triangulate the surface and check if every edge is turnable
    if any(length.(get_facets(polyhedron)).!=3)
        polyhedron = triangulate(polyhedron)
    end

    for edge in get_edges(polyhedron)
        if !isturnable(edge, polyhedron,tol=tol)
            return 0
        end
    end
    return 1
end


"""
poly::Polyhedron
returns a vector of tetrahedra, which union is the Polyhedron poly.
"""
function convexdecomp(poly::Polyhedron; tol::Real=1e-5)::Vector{Polyhedron}
    sol = Polyhedron[]

    triangPoly = triangulate(poly)
    subPoly = deepcopy(triangPoly)

    while !isconvex(subPoly)
        vertexDegrees = map(v -> length(FacesOfVertex(subPoly, v)), 1:length(subPoly.verts))        

        # remove a tetrahedron from subPoly if one exists
        if any(vertexDegrees.==3)
            # println("\nsubPoly:", subPoly)
            # println("vertexDegrees.==3:", [vertexDegrees.==3])
            v = indexin(3, vertexDegrees)[1]
            
            tetraVerts = push!(union(map(e -> setdiff(e, [v]), EdgesOfVertex(subPoly, v))...), v) # determine vertex indices of tetrahedron
            # println("tetraverts:", tetraVerts)
            tetra = Polyhedron(map(w -> subPoly.verts[w], tetraVerts), [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], [[1,2,3], [1,2,4], [1,3,4], [2,3,4]]) # resulting tetrahedron
            # println("found tetra:", tetra)
            push!(sol, tetra)

            # remove v and adjacent edges and facets from subPoly
            subPoly.verts = subPoly.verts[1:end .!= v]
            setdiff!(subPoly.edges, subPoly.edges[map(e -> v in e, subPoly.edges)])
            setdiff!(subPoly.facets, subPoly.facets[map(f -> v in f, subPoly.facets)])

            # add bottom of tetrahedron as facet
            push!(subPoly.facets, tetraVerts[1:3])

            # shift indices of all vertices with index > v
            subPoly.edges = map(e -> map((w -> w>v ? w-1 : w), e), subPoly.edges)
            subPoly.facets = map(f -> map((w -> w>v ? w-1 : w), f), subPoly.facets)

            # println(sol)
            continue
        end

        # edge turn if no tetrahedron was removed
        for edge in subPoly.edges
            # if edge is turnable, cut the resulting tetrahedron out
            if isturnable(edge, subPoly, tol = tol)
                butterfly = subPoly.facets[map(f -> length(Base.intersect(f, edge)) == 2, subPoly.facets)]
                butterflyTips = setdiff(union(butterfly...), edge)

                # println("\nsubpoly: ", subPoly)
                # println("evaluated edge: ", edge)
                # println("corresponding butterfly: ", butterfly)
                setdiff!(subPoly.edges, [edge])
                push!(subPoly.edges, butterflyTips)

                setdiff!(subPoly.facets, butterfly)
                append!(subPoly.facets, [[edge[1], butterflyTips[1], butterflyTips[2]], [edge[2], butterflyTips[1], butterflyTips[2]]])

                tetra = Polyhedron(subPoly.verts[[edge[1], edge[2], butterflyTips[1], butterflyTips[2]]], [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], [[1,2,3], [1,2,4], [1,3,4], [2,3,4]])
                push!(sol, tetra)
                # println("length(sol) = ", length(sol))
                break
            end
        end
    end

    push!(sol, subPoly) # subPoly is now convex

    return sol
end