include("Polyhedron.jl")
include("combinatorics.jl")
include("affine_geometry.jl")
include("polygonal_geometry.jl")


# triangulate surface of Polyhedron
"""
returns a polyhedron containing the vertices and edges of poly such that every facet is triangular.
"""
function triangulate!(poly::Polyhedron; atol = 1e-8)
    coords = get_verts(poly)
    newedges = get_edges(poly)
    newfacets = Vector{Int}[]

    # triangulate every facet of poly by the earcut algorithm
    for facet in get_facets(poly)
        triangulation = earcut3d(coords[facet])
        append!(newfacets, [facet[triang] for triang in triangulation])
        triangulationedges = vcat([ [facet[triang[[1,2]]], facet[triang[[1,3]]], facet[triang[[2,3]]]] for triang in triangulation ]...)
        append!(newedges, triangulationedges)
        newedges = collect.(unique(Set.(newedges)))
    end

    set_edges!(poly, newedges)
    set_facets!(poly, newfacets)
end


function triangulate(poly::Polyhedron; atol = 1e-5)::Polyhedron
    polycopy = deepcopy(poly)
    triangulate!(polycopy, atol = atol)

    return polycopy
end


function outward_normal(poly::Polyhedron{S, T}, facet::Vector{T}) where {S<:Real, T<:Integer}
    @assert facet in get_facets(poly) || reverse(facet) in get_facets(poly) "facet needs to be a facet of poly."
    poly_orient = orient_facets(poly)

    if facet in get_facets(poly_orient)
        f = facet
    else
        f = reverse(facet)
    end

    basis_inds = sort(affinebasis_indices(hcat(get_verts(poly)[f]...)))
    @assert length(basis_inds) == 3

    return normalize(cross(get_verts(poly)[f[basis_inds[2]]] - get_verts(poly)[f[basis_inds[1]]], get_verts(poly)[f[basis_inds[3]]] - get_verts(poly)[f[basis_inds[1]]]))
end


"""
    edgetype(poly::Polyhedron, edge::Vector{<:Int}; atol::Real = 1e-8)

Determine the type of the edge of poly as either "flat", "concave" or "convex".
"""
function edgetype(poly::Polyhedron, edge::Vector{<:Int}; atol::Real = 1e-12)
    @assert edge in get_edges(poly)|| reverse(edge) in get_edges(poly) "edge has to be an edge of poly."
    verts = get_verts(poly)

    poly_orient = orient_facets(poly)
    
    facets = adjfacets(poly_orient, edge)

    function direction(facet, edge)
        # function that determines if edge is oriented forwards (1) or backwards (-1) in facet. 
        inds = indexin(edge, facet)
        if mod1(inds[2] - inds[1], length(facet)) == 1 return 1 end
        return -1
    end

    # edge is forward in one of the facets and backwards in the other as poly_orient is oriented
    inds = indexin(edge, facets[1])
    if direction(facets[1], edge) == 1
        forw = facets[1]
        backw = facets[2]
    else
        forw = facets[2]
        backw = facets[1]
    end

    i = indexin(edge[1], forw)[1]
    j = indexin(edge[1], backw)[1]

    # use determinant to determine if edge forms a right or left system with neighboring edges.
    d = det(hcat(verts[forw[mod1(i-1, length(forw))]] - verts[edge[1]], verts[edge[2]] - verts[edge[1]], verts[backw[mod1(j+1, length(backw))]] - verts[edge[1]]))

    if abs(d) < atol
        return "flat"
    elseif d > 0
        return "convex"
    elseif d < 0
        return "concave"
    end
end


"""
    isflat(poly::Polyhedron, edge::Vector{<:Int})

Determine whether the edge edge of the Polyhedron poly is flat.
"""
function isflatedge(poly::Polyhedron, edge::Vector{<:Int}; atol::Real = 1e-12)
    return edgetype(poly, edge, atol = atol) == "flat"
end


"""
    flattenfacets!(poly::Polyhedron; atol = 1e-5)::Polyhedron

Remove flat edges of the Polyhedron poly.
"""
function flattenfacets!(poly::Polyhedron; atol = 1e-8)
    edgeschecked = []
    
    # @info "edgestocheck: $(edgestocheck)"
    # determine flat edges by calculating the affine dimension of the union of two facets.
    for e in get_edges(poly)
        # @info "e: $(e)"
        if e in edgeschecked
            continue
        end

        @assert length(adjfacets(poly, e)) <= 2 "Edge $(e) is degenerate as it is edge of the facets $(adjfacets(poly, e))"

        if length(adjfacets(poly, e)) < 2
            push!(edgeschecked, e)
            continue
        end

        facet1 = adjfacets(poly, e)[1]
        edges1 = filter(edge -> issubset(edge, facet1), get_edges(poly))
        facet2 = adjfacets(poly, e)[2]
        edges2 = filter(edge -> issubset(edge, facet2), get_edges(poly))
        # @info "facet1: $(facet1)"
        # @info "facet2: $(facet2)"
        # @info "edges1: $(edges1)"
        # @info "edges2: $(edges2)"
        # possible edges in the intersection of the facets. All of those are checked at once since the affine dimension of a facet is 2. 
        # So one of them is flat <=> all of them are flat
        intersect_verts = Base.intersect(facet1, facet2)
        intersect_edges = Base.intersect(edges1, edges2)
        
        append!(edgeschecked, intersect_edges)
        unique!(edgeschecked)

        # if facets span a space of dimension 2 the edges in the intersection can be removed
        # the facets need to be merged
        if isflatedge(poly, e)
            # mark edges in the intersection as checked
            append!(edgeschecked, intersect_edges)


            # remove flat edges from poly
            oldedges = get_edges(poly)
            newedges = setdiff(oldedges, intersect_edges)
            set_edges!(poly, newedges)

            # merge facets
            # verts in the intersection ordered along a path defined by the intersect_edges
            intersect_path = formpath(intersect_verts, intersect_edges)
            # relevant vertices of facet1 and facet2 that are part of the merged facet: the inner vertices of the intersection path are removed
            newfacet1 = setdiff(facet1, intersect_path[2:end-1])
            newfacet2 = setdiff(facet2, intersect_path[2:end-1])
            # edges of facet1 and facet2 that remain after removing flat edges
            relevantedges = symdiff(edges1, edges2)
            mergedfacet = formpath(union(newfacet1, newfacet2), relevantedges)

            oldfacets = get_facets(poly)
            newfacets = setdiff(oldfacets, [facet1, facet2])
            push!(newfacets, mergedfacet)
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
    isturnable(e::Vector{<:Int}, polyhedron::Polyhedron; atol::Real=1e-5)::Bool

Checks whether the edge e is turnable edge of poly. I.e. poly is 
    - a spherical simplicial surface and 
    - the line connecting the wingtips of the butterfly with inner edge e is contained in poly.
Floats with abs value <atol are considered zeros
"""
# TODO: Reihenfolge der Argumente überall glattziehen: immer zuerst Polyhedron übergeben, oder das, worauf sich Methode bezieht?
function isturnable(e::Vector{<:Int}, polyhedron::Polyhedron; atol::Real=1e-5)::Bool
    @assert all(length.(get_facets(polyhedron)).==3) "poly may only have triangle facets"
    @assert in(e, get_edges(polyhedron)) || in(reverse(e), get_edges(polyhedron)) "e has to be an edge of poly"
    
    # adjacent facets to e
    butterfly = adjfacets(polyhedron, e)
    # turned edge
    butterflytips = setdiff(union(butterfly...), e)

    # if there is already an edge connecting the butterflytips, e can't be turned.
    if butterflytips in get_edges(polyhedron) || reverse(butterflytips) in get_edges(polyhedron)
        return 0
    end

    mid = (get_verts(polyhedron)[butterflytips[1]] + get_verts(polyhedron)[butterflytips[2]]) / 2

    # if midpoint of turned edge is inside the polyhedron the edge is turnable
    if inpolyhedron(mid, polyhedron, atol = atol) != 1 # TODO: Testen: macht hier flattenfacets(polyhedron) die Methode wirklich robuster? Wie kann flattenfacets optimiert werden -> sie ist jetzt seeeeehr langsam.
        return 0
    end

    return 1
end


"""
    isconvex(poly::Polyhedron; atol::Real=1e-5)::Bool


determine whether the polyhedron poly is convex. Floats with abs value <atol are considered zeros
"""
function isconvex(poly::Polyhedron; atol::Real=1e-5)::Bool
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
        if !isturnable(edge, polyhedron,atol=atol)
            return 0
        end
    end
    return 1
end



#############################################
# Workaround für tiblöcke!! TODO: debug convex decomposition!
#############################################
function convexdecomp(n1::Integer, n2::Integer, nmerges::Integer)
    @assert mod(n1, nmerges) == 0 "Number of blocks on the side needs to divide n1."
    sol = Polyhedron[]

    if n2 == 3
        sideelem = nprism(3)
        merge!(sideelem, nprism(3), [[3,1,6,4]], [[1,2,4,5]])
        merge!(sideelem, nprism(3), [[2,3,5,6]], [[1,2,4,5]])
    elseif n2 == 4
        sideelem = nprism(4)
        merge!(sideelem, nprism(4), [[2,3,6,7]], [[1,2,5,6]])
        merge!(sideelem, nprism(4), [[4,1,8,5]], [[1,2,5,6]])
    else
        sideelem = nprism(n2)   
    end

    block = nprism(n1)
    push!(sol, block)
    for k in 3:Int(n1 / nmerges):length(get_facets(nprism(n1)))
        preim = get_verts(sideelem)[[n2+1,1,2]]
        push!(preim, preim[1] + cross(preim[2] - preim[1], preim[3] - preim[1]))
        im = get_verts(nprism(n1))[get_facets(nprism(n1))[k][1:3]]
        push!(im, im[1] - cross(im[2] - im[1], im[3] - im[1]))
        
        aff = rigidmap(preim, im)
        newsideelem = deepcopy(sideelem)
        set_verts!(newsideelem, aff.(get_verts(newsideelem)))
        push!(sol, newsideelem)
    end

    return sol
end

"""
poly::Polyhedron
returns a vector of tetrahedra, which union is the Polyhedron poly.
"""
function convexdecomp(poly::Polyhedron; atol::Real=1e-5)::Vector{Polyhedron}
    #############################################
    # Workaround für tiblöcke!! TODO: debug convex decomposition!
    #############################################
    for n1 in 3:20
        for n2 in 3:10
            for nmerges in (1:n1)[mod.(n1, 1:n1) .== 0]
                if iscongruent(tiblock(n1, n2, nmerges), poly)
                    sol = Polyhedron[]

                    if n2 == 3
                        sideelem = nprism(3)
                        merge!(sideelem, nprism(3), [[3,1,6,4]], [[1,2,4,5]])
                        merge!(sideelem, nprism(3), [[2,3,5,6]], [[1,2,4,5]])
                    elseif n2 == 4
                        sideelem = nprism(4)
                        merge!(sideelem, nprism(4), [[2,3,6,7]], [[1,2,5,6]])
                        merge!(sideelem, nprism(4), [[4,1,8,5]], [[1,2,5,6]])
                    else
                        sideelem = nprism(n2)   
                    end

                    block = nprism(n1)
                    for k in 3:Int(n1 / nmerges):length(get_facets(nprism(n1)))
                        preim = get_verts(sideelem)[[n2+1,1,2]]
                        push!(preim, preim[1] + cross(preim[2] - preim[1], preim[3] - preim[1]))
                        im = get_verts(nprism(n1))[get_facets(nprism(n1))[k][1:3]]
                        push!(im, im[1] - cross(im[2] - im[1], im[3] - im[1]))
                        
                        aff = rigidmap(preim, im)
                        newsideelem = deepcopy(sideelem)
                        set_verts!(newsideelem, aff.(get_verts(newsideelem)))
                        push!(sol, newsideelem)
                    end
                    push!(sol, block)

                    rigid = rigidmap(get_verts(tiblock(n1, n2, nmerges)), get_verts(poly))
                    for component in sol
                        set_verts!(component, rigid.(get_verts(component)))
                    end

                    return sol
                end
            end
        end
    end


    sol = Polyhedron[]

    triangPoly = triangulate(poly)
    subPoly = deepcopy(triangPoly)

    display(plot(subPoly, labels = true, width = 400, height = 300))
    if any(length.(get_facets(subPoly)) .> 3)
        @warn "facets with too many vertices: $(get_facets(subPoly)[length.(get_facets(subPoly)) .> 3])"
    end
    while !isconvex(subPoly)
        if any(length.(get_facets(subPoly)) .> 3)
            @warn "facets with too many vertices: $(get_facets(subPoly)[length.(get_facets(subPoly)) .> 3])"
        end
        vertexDegrees = map(v -> length(incfacets(subPoly, v)), 1:length(subPoly.verts))        

        # remove a tetrahedron from subPoly if one exists
        if any(vertexDegrees.==3)
            # println("\nsubPoly:", subPoly)
            # println("vertexDegrees.==3:", [vertexDegrees.==3])
            v = indexin(3, vertexDegrees)[1]
            @info "v: $(v)"
            @info "faces of v: $(incfacets(subPoly, v))"
            
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
            @info "bottom tetrahedron: $(tetraVerts[1:3])"

            # shift indices of all vertices with index > v
            subPoly.edges = map(e -> map((w -> w>v ? w-1 : w), e), subPoly.edges)
            subPoly.facets = map(f -> map((w -> w>v ? w-1 : w), f), subPoly.facets)

            # println(sol)
            display(plot(subPoly, labels = true, width = 400, height = 300))
            continue
        end

        # edge turn if no tetrahedron was removed
        for edge in subPoly.edges
            # if edge is turnable, cut the resulting tetrahedron out
            if isturnable(edge, subPoly, atol = atol) && !isflatedge(subPoly, edge)
                @info "edge: $(edge)"
                butterfly = filter(f -> length(Base.intersect(f, edge)) == 2, get_facets(subPoly))
                @info "butterfly: $(butterfly)"
                butterflyTips = setdiff(union(butterfly...), edge)
                @info "butterfly tips: $(butterflyTips)"

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
                display(plot(subPoly, labels = true, width = 400, height = 300))
                break
            end
        end
    end

    push!(sol, subPoly) # subPoly is now convex

    return sol
end


