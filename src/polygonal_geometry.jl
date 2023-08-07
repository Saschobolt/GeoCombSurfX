include("affine_geometry.jl")

function intriang3d(triang::Vector{<:Vector{<:Real}}, p::Vector{<:Real})
    # https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html
    @assert all(length.(triang) .== 3) "triang has to consist of 3-vectors."
    @assert length(triang) == 3 "triang needs to be a list of 3 points."
    @assert length(p) == 3 "p has to be a point in 3-space."
    @assert affinedim(triang) == 2 "triang is degenerate."

    # determine whether p lies on the boundary of triang
    for i in 1:3
        for j in (i+1):3
            if affinedim([triang[i], triang[j], p]) == 1 && dot(triang[i] - p, triang[j] - p) < 0
                return -1
            end
        end
    end

    # translate p into origin
    a = triang[1] - p
    b = triang[2] - p
    c = triang[3] - p

    # compute normals of triangles pab, pbc, pca
    u = cross(a,b)
    v = cross(b,c)
    w = cross(c,a)

    # determine whether the normals face the same direction.
    # then the triangles are wound the same and thus p is inside the triangle
    if dot(u,v) < 0.0
        return 0
    end

    if dot(u,w) < 0.0
        return 0
    end

    return 1
end


function earcut3d(polygon::Vector{<:Vector{<:Real}}, atol = 1e-8)
    # earcut algorithm: https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
    # https://www.mathematik.uni-marburg.de/~thormae/lectures/graphics1/code/JsCoarseImg/EarCutting.html
    @assert all(length.(polygon) .== 3) "polygon needs to be 3d Polygon."
    # TODO: assert that polygon is simple

    if polygon[1] == polygon[end]
        coords = polygon[1:end-1]
    else
        coords = polygon
    end

    remaining_verts = collect(1:length(polygon))
    sol = []

    # @info "remaining_verts: $(remaining_verts)"

    # return the neighbors of vertex v in the polygon indexed by subfacet
    function neighbors(v::Int)
        i = indexin(v, remaining_verts)[1]
        neighbor1 = remaining_verts[mod1(i-1, length(remaining_verts))]
        neighbor2 = remaining_verts[mod1(i+1, length(remaining_verts))]
        return neighbor1, neighbor2
    end

    # normal vector of polygon Plane
    n = normalvec(coords[remaining_verts])

    # signed angles at vertices
    angles  = [signedangle3d(coords[neighbors(v)[1]] - coords[v], coords[v] - coords[neighbors(v)[2]], n, atol = 1e-12) for v in remaining_verts]
    angles[abs.(angles) .< 1e-12] .= 0

    # sign of angles at convex vertices is the same as the sum of angle array (2pi or -2pi)
    convexsign = sign(sum(angles))

    # boolean to determine whether a vertex is convex
    function isconvex(v::Int)
        neighbor1 = neighbors(v)[1]
        neighbor2 = neighbors(v)[2]
        angle = signedangle3d(coords[neighbor1] - coords[v], coords[v] - coords[neighbor2], n, atol = 1e-12)
        # @info "check, if $(v) is convex; neighbor1: $(neighbor1), neighbor2: $(neighbor2), angle: $(angle)"
        if abs(angle) < 1e-12
            return false
        end
        return sign(angle) == convexsign
    end

    # boolean to determine whether a vertex is reflex
    function isreflex(v::Int)
        return !isconvex(v)
    end

    # initialize vertex arrays
    convex  = remaining_verts[isconvex.(remaining_verts)]
    reflex  = remaining_verts[isreflex.(remaining_verts)]

    # boolean to check whether convex vertex v is an ear, i.e.
    # - all reflex vertices lie outside the triangle spanned by v and its neighbors 
    # - the affine dimension of the polygon after removing v is still 2
    function isear(v::Int)
        neighbor1, neighbor2 = neighbors(v)
        triangle = [coords[neighbor1], coords[v], coords[neighbor2]]

        # @info "check, if $(v) is an ear; neighbor1: $(neighbor1), neighbor2: $(neighbor2)"

        return all([intriang3d(triangle, coords[w]) == 0 for w in setdiff(reflex, [neighbor1, neighbor2])]) && affinedim(coords[setdiff(remaining_verts, [v])]) == 2
    end

    # calculate the ration of the area and the circumference of an ear: area / circumference.
    function acratio(v::Int)
        neighbor1, neighbor2 = neighbors(v)
        triangle = [coords[neighbor1], coords[neighbor2], coords[v]]
        side1 = triangle[2] - triangle[1]
        side2 = triangle[3] - triangle[2]
        side3 = triangle[1] - triangle[3]

        # volume and circumference of the triangle
        volume = det(hcat(side1, -side3, normalize(cross(side1, -side3)))) / 2
        circumference = sum(norm.([side1, side2, side3]))

        return volume / circumference
    end

    ears = sort(filter(v -> isear(v), convex), by = acratio)

    # @info "convex: $(convex)"
    # @info "reflex: $(reflex)"
    # @info "ears: $(ears)"

    while length(remaining_verts) > 3
        # remove ears one at a time
        v = pop!(ears)
        # @info "removing ear $(v)"
        neighbor1 = neighbors(v)[1]
        neighbor2 = neighbors(v)[2]
        ear = [neighbor1, v, neighbor2]

        # ear is a triangle in the triangulation
        push!(sol, ear)

        # remove v from subfacet and vertex arrays
        setdiff!(remaining_verts, [v])
        setdiff!(convex, [v])


        # reevaluate neighbors: convex stays convex, ears don't have to stay ears, reflex can become convex or ear
        for w in [neighbor1, neighbor2]
            # @info "investigating neighbor $(w)"
            if w in reflex
                if isconvex(w)
                    setdiff!(reflex, [w])
                    push!(convex, w)
                    isear(w) ? push!(ears, w) : continue
                end
            end

            if w in convex && !(w in ears)
                if isear(w) 
                    push!(ears, w)
                end

                continue
            end

            if w in ears
                if !isear(w)
                    setdiff!(ears, [w])
                end

                continue
            end
        end

        sort!(ears, by = acratio)

        # @info "remaining verts: $(remaining_verts)"
        # @info "new convex: $(convex)"
        # @info "new reflex: $(reflex)"
        # @info "new ears: $(ears)"
    end

    # if the subfacet is a triangle, add it to the facet list and continue with the next facet
    push!(sol, remaining_verts)

    return sol
end



"""
    inpolygon3d(poly::Vector{<:Vector{<:Real}}, p::Vector{<:Real}; atol = 1e-8)

Determine whether the point p lies inside the polygon poly.

1: p is inside polygon
0: p is outside polygon
-1: p lies on the boundary of polygon
"""
function inpolygon3d(poly::Vector{<:Vector{<:Real}}, p::Vector{<:Real}; atol = 1e-8)
    @assert all(length.(poly) .== 3) "triang has to consist of 3-vectors."
    @assert length(poly) == 3 "triang needs to be a list of 3 points."
    @assert length(p) == 3 "p has to be a point in 3-space."
    @assert affinedim(poly) == 2 "triang is degenerate."

    poly_triang = [poly[triangle] for triangle in earcut3d(poly, atol = atol)]
    if any([intriang3d(triang, p) == 1 for triang in poly_triang])
        return 1
    end
    # TODO: return -1 if p lies on the boundary of the polygon
    return 0
end



# """
# returns 1 if the 3d point p lies in the polygon poly embedded into R^3, -1 if it lies on its boundary and 0 otherwise
# """
# function inpolygon3d(p::AbstractVector, poly::AbstractVector; tol::Real=1e-5)
#     d = 3
#     @assert length(p) == 3 "Only 3d case."
#     @assert all([length(v) == 3 for v in poly]) "Only 3d case."

#     E = Plane(poly, tol = tol)
#     point = deepcopy(p)

#     preim = [E.point, E.point + E.vectors[1], E.point + E.vectors[2], E.point + normalize(cross(E.vectors[1], E.vectors[2]))]
#     im = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]

#     # affine map so that the affine basis of the polygon is mapped to the xy plane
#     aff = affinemap(preim, im, atol = tol)

#     # transformed polygon and point
#     polytransformed = aff.(poly)
#     pointtransformed = aff(point)
    

#     if abs(pointtransformed[3]) < tol
#          # inpolygon projects onto the xy plane and decides the problem.
#          # If p and poly lie in the same plane, then this is enough to decide the problem
#         return inpolygon(pointtransformed, polytransformed)
#     else
#         # p and poly don't lie in the same plane
#         return 0 
#     end
# end