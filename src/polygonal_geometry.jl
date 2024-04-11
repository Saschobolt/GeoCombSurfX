include("affine_geometry.jl")

"""
    to_xyplane_map(polygon::AbstractMatrix{T}; atol::Real = 1e-8) where T<:Real

Map that transforms the polygon so that it lies in the xy plane. Triangles will be transformed to the upper halfspace. 
Linear part of map has determinant 1 and thus retains orientation. 
"""
function to_xyplane_map(polygon::AbstractMatrix{T}; atol::Real=1e-8) where {T<:Real}
    # https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane
    n = size(polygon)[2]

    if size(polygon)[1] == 2
        return vcat(polygon, repeat([0], 1, n))
    end

    @assert affinedim(polygon, atol=atol) == 2 "Polygon doesn't span a plane."
    @assert size(polygon)[1] == 3 "only implemented for 3d polygons."

    # random translation vector to apply to the polygon so that no vertex lies at the origin
    r = rand(promote_type(Float64, T), size(polygon)[1])
    polygon_trafo = copy(polygon) + repeat(r, 1, n)

    j = findfirst(i -> affinedim(polygon_trafo[:, [i, mod1(i + 1, n), mod1(i + 2, n)]], atol=atol) == 2, collect(1:n))
    s = -polygon_trafo[:, j]
    polygon_trafo = copy(polygon_trafo) + repeat(s, 1, n)
    triang = polygon_trafo[:, [j, mod1(j + 1, n), mod1(j + 2, n)]]
    u = normalize(triang[:, 2])
    w = normalize(cross(u, triang[:, 3]))
    v = cross(u, w)
    M = transpose(hcat(u, v, w))
    M = det(M) * M

    f = x::AbstractVecOrMat -> typeof(x) <: AbstractVector ? M * (x + r + s) : M * (x + repeat(r, 1, n) + repeat(s, 1, n))
    poly_trafo = f(polygon)

    # want that triangle lies in first or second quadrant (first vertex at [0,0,0], second at [d,0,0], third at [e,f,0]). Calculate rotation matrix that achieves that.
    rot = [1 0 0; 0 1 0; 0 0 1]
    if poly_trafo[:, 2][1] < 0
        rot = [-1 0 0; 0 -1 0; 0 0 1]
        poly_trafo = rot * poly_trafo
    end
    if poly_trafo[:, 3][2] < 0
        rot = [1 0 0; 0 -1 0; 0 0 -1] * rot
    end

    f = x::AbstractVecOrMat -> typeof(x) <: AbstractVector ? rot * M * (x + r + s) : rot * M * (x + repeat(r, 1, n) + repeat(s, 1, n))

    return f
end

"""
    to_xyplane(poly::AbstractMatrix{<:Real})

Transform the polygon so that it lies in the xy plane. Triangles will be transformed to upper halfspace. Orientation is retained.
"""
function to_xyplane(polygon::AbstractMatrix{<:Real}; atol::Real=1e-8)
    # transform polygon so that it lies in the xy plane
    return to_xyplane_map(polygon, atol=atol)(polygon)
end

"""
    intriang(triang::AbstractMatrix{<:Real}, p::Vector{<:Real})

Determine whether the point p lies inside the triangle triang.
"""
function intriang(triang::AbstractMatrix{<:Real}, p::AbstractVector{<:Real}; atol=1e-8)
    if affinedim(triang; atol=atol) < 2
        display(triang)
        display(affinedim(triang))
        error("triang is degenerate.")
    end

    if size(triang)[1] == 3
        # if point and triangle don't lie in the same plane, p is not inside triang
        if affinedim(hcat(triang, p); atol=atol) > 2
            return 0
        end

        # transform triang and p so that triang lies in the xy plane
        f = to_xyplane_map(triang, atol=atol)
        triang_trafo = f(triang)[1:2, :]
        p_trafo = f(p)[1:2]
        return intriang(triang_trafo, p_trafo; atol=atol)
    end

    # determine baricentric coordinates of p with respect to the triangle
    leg1 = triang[:, 2] - triang[:, 1]
    leg2 = triang[:, 3] - triang[:, 1]

    A = hcat(leg1, leg2)
    b = p - triang[:, 1]
    s = A \ b

    if all(s .> atol) && abs(sum(s) - 1) > atol
        return 1
    elseif any(abs.(s) .< atol) || (all(s .> atol) && abs(sum(s) - 1) < atol)
        return -1
    else
        return 0
    end
end

function intriang(triang::Vector{<:Vector{<:Real}}, p::Vector{<:Real}; atol=1e-8)
    return intriang(hcat(triang...), p; atol=atol)
end


"""
    is_ccw(polygon::Matrix{<:Real}, n::Vector{<:Real}; atol = 1e-12)

Determine whether the orientation of the polygon is counterclockwise. Clockwise means, that the vertices follow a negative rotation around the normal vector n.
"""
function is_ccw(polygon::AbstractMatrix{<:Real}, n::Vector{<:Real}; atol=1e-12)
    # https://math.stackexchange.com/questions/2152623/determine-the-order-of-a-3d-polygon
    if polygon[:, 1] == polygon[:, end]
        coords = polygon[:, 1:end-1]
    else
        coords = polygon
    end

    m = size(coords)[2]

    @assert all([dot(n, coords[:, mod1(i + 1, m)] - coords[:, mod1(i, m)]) == 0 for i in 1:m]) "n is not normal to the polygon."

    s = sum([cross(coords[:, mod1(i, m)], coords[:, mod1(i + 1, m)]) for i in 1:m])

    if dot(s, n) > 0
        return true
    else
        return false
    end
end


# function earcut3d(polygon::AbstractMatrix{<:Real}; atol = 1e-8)
#     # earcut algorithm: https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
#     # https://www.mathematik.uni-marburg.de/~thormae/lectures/graphics1/code/JsCoarseImg/EarCutting.html
#     @assert size(polygon)[1] == 3 "polygon needs to be 3d Polygon."
#     # TODO: assert that polygon is simple

#     if polygon[:,1] == polygon[:,end]
#         coords = polygon[:,1:end-1]
#     else
#         coords = polygon
#     end

#     remaining_verts = collect(1:size(coords)[2])
#     if length(remaining_verts) == 3
#         return [remaining_verts]
#     end
#     sol = Vector{Int}[]

#     # @info "remaining_verts: $(remaining_verts)"

# # return the neighbors of vertex v in the polygon indexed by subfacet
# function neighbors(v::Int)
#     i = indexin(v, remaining_verts)[1]
#     neighbor1 = remaining_verts[mod1(i-1, length(remaining_verts))]
#     neighbor2 = remaining_verts[mod1(i+1, length(remaining_verts))]
#     return neighbor1, neighbor2
# end

#     # normal vector of polygon Plane
#     n = normalvec(coords[:,remaining_verts])

#     # signed angles at vertices
#     angles  = [signedangle3d_right(coords[:,neighbors(v)[1]] - coords[:,v], coords[:,v] - coords[:,neighbors(v)[2]], n, atol = atol) for v in remaining_verts]
#     angles[abs.(angles) .< atol] .= 0

#     # sign of angles at convex vertices is the same as the sum of angle array (2pi or -2pi)
#     convexsign = sign(sum(angles))

#     # boolean to determine whether a vertex is convex
#     function isconvex(v::Int)
#         neighbor1 = neighbors(v)[1]
#         neighbor2 = neighbors(v)[2]

#         if affinedim(coords[:, [neighbor1, v, neighbor2]]) < 2
#             return false
#         end

#         n1 = normalvec(coords[:,[neighbor1, v, neighbor2]])
#         if dot(n1, n) < 0
#             n1 = -n1
#         end

#         angle = signedangle3d_right(coords[:,neighbor1] - coords[:,v], coords[:,v] - coords[:,neighbor2], n, atol = atol)
#         # @info "check, if $(v) is convex; neighbor1: $(neighbor1), neighbor2: $(neighbor2), angle: $(angle)"
#         if abs(angle) < atol
#             return false
#         end
#         return sign(angle) == convexsign
#     end

#     # boolean to determine whether a vertex is reflex
#     function isreflex(v::Int)
#         return !isconvex(v)
#     end

#     # initialize vertex arrays
#     convex  = remaining_verts[isconvex.(remaining_verts)]
#     reflex  = remaining_verts[isreflex.(remaining_verts)]

#     # boolean to check whether convex vertex v is an ear, i.e.
#     # - all reflex vertices lie outside the triangle spanned by v and its neighbors 
#     # - the affine dimension of the polygon after removing v is still 2
#     function isear(v::Int)
#         neighbor1, neighbor2 = neighbors(v)
#         triangle = coords[:, [neighbor1, v, neighbor2]]

#         # @info "check, if $(v) is an ear; neighbor1: $(neighbor1), neighbor2: $(neighbor2)"

#         # if the triangle is degenerate, v can't be an ear
#         if affinedim(triangle, atol = atol) < 2
#             return false
#         end

#         return all([intriang3d(triangle, coords[:,w]) == 0 for w in setdiff(reflex, [neighbor1, neighbor2])]) && affinedim(coords[:, setdiff(remaining_verts, [v])]) == 2
#     end

#     # calculate the ration of the area and the circumference of an ear: area / circumference.
#     function acratio(v::Int)
#         neighbor1, neighbor2 = neighbors(v)
#         triangle = [coords[:,neighbor1], coords[:,neighbor2], coords[:,v]]
#         side1 = triangle[2] - triangle[1]
#         side2 = triangle[3] - triangle[2]
#         side3 = triangle[1] - triangle[3]

#         # volume and circumference of the triangle
#         volume = det(hcat(side1, -side3, normalize(cross(side1, -side3)))) / 2
#         circumference = sum(norm.([side1, side2, side3]))

#         return volume / circumference
#     end

#     ears = sort(filter(v -> isear(v), convex), by = acratio)

#     # @info "convex: $(convex)"
#     # @info "reflex: $(reflex)"
#     # @info "ears: $(ears)"

#     while length(remaining_verts) > 3
#         # remove ears one at a time
#         v = pop!(ears)
#         # @info "removing ear $(v)"
#         neighbor1 = neighbors(v)[1]
#         neighbor2 = neighbors(v)[2]
#         ear = [neighbor1, v, neighbor2]

#         # ear is a triangle in the triangulation
#         push!(sol, ear)

#         # remove v from subfacet and vertex arrays
#         setdiff!(remaining_verts, [v])
#         setdiff!(convex, [v])


#         # reevaluate neighbors: convex stays convex, ears don't have to stay ears, reflex can become convex or ear
#         for w in [neighbor1, neighbor2]
#             # @info "investigating neighbor $(w)"
#             if w in reflex
#                 if isconvex(w)
#                     setdiff!(reflex, [w])
#                     push!(convex, w)
#                     isear(w) ? push!(ears, w) : continue
#                 end
#             end

#             if w in convex && !(w in ears)
#                 if isear(w) 
#                     push!(ears, w)
#                 end

#                 continue
#             end

#             if w in ears
#                 if !isear(w)
#                     setdiff!(ears, [w])
#                 end

#                 continue
#             end
#         end

#         sort!(ears, by = acratio)

#         # @info "remaining verts: $(remaining_verts)"
#         # @info "new convex: $(convex)"
#         # @info "new reflex: $(reflex)"
#         # @info "new ears: $(ears)"
#     end

#     # if the subfacet is a triangle, add it to the facet list and continue with the next facet
#     push!(sol, remaining_verts)

#     return sol
# end

function earcut3d(polygon::AbstractMatrix{<:Real}; atol=1e-8)
    # transform polygon so that it lies in the xy plane
    # https://math.stackexchange.com/questions/856666/how-can-i-transform-a-3d-triangle-to-xy-plane
    polygon_trafo = to_xyplane(polygon, atol=atol)[1:2, :]
    if polygon_trafo[:, end] == polygon_trafo[:, 1]
        polygon_trafo = polygon_trafo[:, 1:end-1]
    end

    # return the neighbors of vertex v in the polygon indexed by subfacet
    # TODO: Use a graph structure to store the polygon. Graphs.jl has efficient implementation for neighbor calculations.
    function neighbors(v::Int)
        i = indexin(v, remaining)[1]
        neighbor1 = remaining[mod1(i - 1, length(remaining))]
        neighbor2 = remaining[mod1(i + 1, length(remaining))]
        return neighbor1, neighbor2
    end

    n = size(polygon_trafo)[2]
    remaining = collect(1:n)
    angles = [signedangle2d(polygon_trafo[:, i] - polygon_trafo[:, neighbors(i)[1]], polygon_trafo[:, neighbors(i)[2]] - polygon_trafo[:, i]) for i in remaining]
    # @info "angles: $(angles)"
    sum_signed_angles = sum(angles)
    # @info "sum_signed_angles: $(sum_signed_angles)"
    convex_sign = sign(sum_signed_angles)

    function isconvex(v)
        n1, n2 = neighbors(v)
        angle = signedangle2d(polygon_trafo[:, v] - polygon_trafo[:, n1], polygon_trafo[:, n2] - polygon_trafo[:, v])

        # flat edges are not convex
        if abs(angle) < atol || abs(abs(angle) - pi) < atol
            return false
        end
        return sign(angle) == convex_sign
    end

    convex = filter(v -> isconvex(v), remaining)

    function isreflex(v)
        return !isconvex(v)
    end

    reflex = setdiff(remaining, convex)

    function isear(v)
        # v is an ear, if no reflex vertex lies in the triangle spanned by v and its neighbors in the subfacet
        n1, n2 = neighbors(v)

        return all([intriang(polygon_trafo[:, [n1, v, n2]], polygon_trafo[:, w]; atol=atol) == 0 for w in setdiff(reflex, [n1, n2])])
    end

    ears = filter(v -> isear(v), convex)

    # ratio of area to circumference of the triangle consisting of an ear and its neighbors
    function acratio(v)
        n1, n2 = neighbors(v)
        area = det(hcat(polygon_trafo[:, n1] - polygon_trafo[:, v], polygon_trafo[:, n2] - polygon_trafo[:, v])) / 2
        circ = sum([dist(polygon_trafo[:, n1], polygon_trafo[:, v]), dist(polygon_trafo[:, n2], polygon_trafo[:, v]), dist(polygon_trafo[:, n1], polygon_trafo[:, n2])])
        return area / circ
    end

    # sort ears by acratio so that the first entry in ears is the ear with biggest acratio
    sort!(ears, by=acratio)

    sol = Vector{Int}[]

    # iteratively clip ears of subfacet
    while length(remaining) > 3
        # @info "remaining: $(remaining)"
        # @info "reflex: $(reflex)"
        # @info "convex: $(convex)"
        # @info "ears: $(ears)"

        v = pop!(ears)
        n1, n2 = neighbors(v)

        setdiff!(remaining, [v])
        setdiff!(convex, [v])

        push!(sol, [n1, v, n2])

        # reevaluate neighbors: convex stays convex. ear doesn't have to stay ear. reflex can become convex and ear
        for w in [n1, n2]
            if w in reflex
                if isconvex(w)
                    setdiff!(reflex, [w])
                    push!(convex, w)
                    isear(w) ? push!(ears, w) : continue
                end
            end
        end

        for w in [n1, n2]
            if w in convex
                if w in ears
                    isear(w) ? continue : setdiff!(ears, w)
                else
                    isear(w) ? push!(ears, w) : continue
                    sort!(ears, by=acratio)
                end
            end
        end
    end

    push!(sol, remaining)
    return sol
end



"""
    inpolygon3d(poly::AbstractMatrix{<:Real}, p::Vector{<:Real}; atol = 1e-8)

Determine whether the point p lies inside the polygon poly.

1: p is inside polygon
0: p is outside polygon
-1: p lies on the boundary of polygon
"""
function inpolygon3d(poly::AbstractMatrix{<:Real}, p::Vector{<:Real}; atol=1e-8)
    @assert size(poly)[1] == 3 "poly has to be a 3×n-matrix, but is of size $(size(poly))."
    if poly[:, end] == poly[:, 1]
        coords = poly[:, 1:end-1]
    else
        coords = poly
    end
    n = size(coords)[2]

    for i in 1:size(poly)[2]
        if (max(abs.(p - coords[:, i])...) < atol) || (affinedim([coords[:, i], coords[:, mod1(i + 1, n)], p]; atol=atol) == 1 && dot(coords[:, i] - p, coords[:, mod1(i + 1, n)] - p) <= 0)
            return -1
        end
    end

    poly_triang = [poly[:, triangle] for triangle in earcut3d(poly, atol=atol)]
    for triang in poly_triang
        if intriang(triang, p; atol=atol) != 0
            return 1
        end
    end

    return 0
end