"""
    merge!(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}})

Merge the two polyhedra poly1 and poly2 such that the corresponding vertices of the faces in faces1 and faces2 align. 
The aligned faces are deleted in the resulting polyhedron so that the result itself is a 2-dimensional piecewise linear manifold.
Real values < atol are considered zero.
"""
function merge!(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}}; atol::Real=1e-8)
    @assert all([Set(f) in Set.(get_facets(poly1)) for f in facets1]) "faces1 needs to consist of facets of poly1."
    @assert all([Set(f) in Set.(get_facets(poly2)) for f in facets2]) "faces2 needs to consist of facets of poly2."

    @assert dimension(poly1) == dimension(poly2) "poly1 and poly2 need to be embedded into the same space."

    # matrix with columns describing which vertex indices will be identified in the result
    ind_to_align = unique(hcat(vcat(facets1...), vcat(facets2...)), dims=1)

    @assert length(unique(ind_to_align[:, 1])) == length(ind_to_align[:, 1]) "Merging isn't well defined"
    @assert length(unique(ind_to_align[:, 2])) == length(ind_to_align[:, 2]) "Merging isn't well defined"

    # vertices which are aligned
    verts_to_align1 = get_verts(poly1)[:, ind_to_align[:, 1]]
    verts_to_align2 = get_verts(poly2)[:, ind_to_align[:, 2]]

    for i in 1:length(verts_to_align1)
        for j in i+1:size(verts_to_align2)[2]
            if abs(dist(verts_to_align1[:, i], verts_to_align1[:, j]) - dist(verts_to_align2[:, i], verts_to_align2[:, j])) >= atol
                error("Polyhedra cannot be merged. Distance between Vertex $(ind_to_align[i,1]) and $(ind_to_align[j,1]) of poly1 is $(dist(verts_to_align1[:,i], verts_to_align1[:,j])), but the distance between vertex $(ind_to_align[i,2]) and $(ind_to_align[j,2]) of poly2 is $(dist(verts_to_align2[:, i], verts_to_align2[:, j]))")
            end
        end
    end
    # affine preimage basis and image basis of rigid mapping merging the polyhedra 
    affbasis_indices = affinebasis_indices(verts_to_align2, atol=atol)
    preim = verts_to_align2[:, affbasis_indices]
    im = verts_to_align1[:, affbasis_indices]
    @assert size(preim)[2] >= 3 "Polyhedra cannot be merged along degenerate faces."

    # if facets span a space of dimension 2 add the normal vector of the planes to obtain a basis
    if size(preim)[2] == 3
        preim = hcat(preim, preim[:, 1] + cross(preim[:, 2] - preim[:, 1], preim[:, 3] - preim[:, 1]))
        im = hcat(im, im[:, 1] - cross(im[:, 2] - im[:, 1], im[:, 3] - im[:, 1]))
    end

    # rigid map mapping facets2 to facets1
    tau = rigidmap(preim, im, atol=atol)

    # function that maps the index of a vertex of poly2 to the resulting vertex after merging
    function index_map2(i::Integer)
        if i in ind_to_align[:, 2]
            return ind_to_align[indexin(i, ind_to_align[:, 2])[1], 1]
        end

        sol = i + size(get_verts(poly1))[2] - sum(map(j -> i > j ? 1 : 0, ind_to_align[:, 2]))
        return sol
    end

    # vertices of resulting polyhedron
    sol_verts1 = get_verts(poly1)
    sol_verts2 = tau(get_verts(poly2))[:, map(index -> !(index in ind_to_align[:, 2]), 1:size(get_verts(poly2))[2])]
    sol_verts = hcat(sol_verts1, sol_verts2)

    # edges of resulting polyhedron
    sol_edges1 = get_edges(poly1)
    sol_edges2 = map(e -> index_map2.(e), get_edges(poly2))

    sol_edges = deepcopy(sol_edges1)
    for e in sol_edges2
        if filter(ee -> Set(ee) == Set(e), sol_edges) == []
            push!(sol_edges, e)
        end
    end

    # facets of resulting polyhedron
    sol_facets1 = get_facets(poly1)[map(f -> !(Set(f) in Set.(facets1)), get_facets(poly1))]
    sol_facets2 = get_facets(poly2)[map(f -> !(Set(f) in Set.(facets2)), get_facets(poly2))]
    sol_facets2 = [index_map2.(f) for f in sol_facets2]
    sol_facets = vcat(sol_facets1, sol_facets2)

    # update attributes of poly1
    poly1.verts = sol_verts
    poly1.edges = sol_edges
    poly1.facets = sol_facets
    set_halfedges!(poly1)
end


"""
    merge(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}})

Merge the two polyhedra poly1 and poly2 such that the corresponding vertices of the faces in faces1 and faces2 align. 
The aligned faces are deleted in the resulting polyhedron so that the result itself is a 2-dimensional piecewise linear manifold.
Real values < atol are considered zero.
"""
function merge(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}}; atol::Real=1e-8)
    poly1_copy = deepcopy(poly1)
    merge!(poly1_copy, poly2, facets1, facets2, atol=atol)

    return poly1_copy
end
