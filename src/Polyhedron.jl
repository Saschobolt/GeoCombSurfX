using LinearAlgebra
using PolygonOps

mutable struct Polyhedron
    verts::Vector{Vector{T}} where T<:Number # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{Vector{Int}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{Vector{Int}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
end

struct Ray
    point::Vector{T} where T<:Number
    vector::Vector{T} where T<:Number
end

struct Plane
    point::Vector{T} where T<:Number
    vectors::Vector{Vector{T}} where T<:Number
end


"""
A::Matrix
cond::Float64
Returns an orthonormal basis for the columnspace of the matrix A using svd. Singular values with abs < cond are treated as 0.
"""
function colspace(A::Matrix{T}, cond::Float64 = 1e-5)::Vector{Vector{Float64}} where T<:Number
    F = svd(A)
    return [c[:] for c in eachcol(F.U[:, findall(>(cond), abs.(F.S))])]
end


"""
A::Matrix
cond::Float64
Returns an orthonormal basis for the nullspace of the matrix A using svd. Singular values with abs < cond are treated as 0.
"""
function nullspace(A::Matrix{T}, cond::Float64 = 1e-5)::Vector{Vector{Float64}} where T<:Number
    F = svd(A)
    return [c[:] for c in eachcol(F.V[:, findall(<(cond), abs.(F.S))])]
end

"""
returns the plane, in which polygon lies.
"""
function plane(polygon::Vector{Vector{T}})::Plane where T<:Number
    @assert polygon[1] == polygon[end] "First and last vertex of polygon have to equal."

    v = polygon[1]
    A = hcat(map(w -> w-v, polygon[2:end-1])...)
    return Plane(v, colspace(A))
end

"""
returns the intersection between a ray and a plane if it exists.
"""
function intersect(ray::Ray, plane::Plane)
    vr = ray.point
    vp = plane.point

    A = hcat(union(plane.vectors, [-ray.vector])...)
    if length(colspace(A)) < 3
        throw(ErrorException("ray and plane are parallel."))
    end
    
    coeffs = A\(vr - vp)
    if coeffs[3] < 0
        throw(ErrorException("ray and plane don't intersect"))
    end

    return vr + coeffs[3] * ray.vector
end


"""
returns 1 if the 3d point p lies in the polygon poly embedded into R^3, -1 if it lies on its boundary and 0 otherwise
"""
function inpolygon3d(p::AbstractVector, poly::AbstractVector, tol::Float64=1e-5)
    E = plane(poly)
    point = deepcopy(p)

    trafo = hcat(union(E.vectors, [normalize(cross(E.vectors...))])...)' # transformation so that E is parallel to the xy plane
    E = Plane(trafo * E.point, map(v -> trafo * v, E.vectors))
    point = trafo * point
    polytransformed = map(v -> trafo * v, poly)

    # display(hcat(E.vectors...))
    # println((hcat(E.vectors...) * (hcat(E.vectors...)\(point-E.point))) - (point-E.point))

    if max(abs.((hcat(E.vectors...) * (hcat(E.vectors...)\(point-E.point))) - (point-E.point))...) < tol
        return inpolygon(point, polytransformed) # inpolygon projects onto the xy plane and decides the problem. If p and poly lie in the same plane, then this is enough to decide the problem
    else
        return 0 # p and poly don't lie in the same plane
    end
end


# decide whether point is inside of polyhedron
"""
    Randomized algorithm to check whether a point is contained in a polyhedron.
"""
function inpolyhedron(point::Vector{T}, poly::Polyhedron, tol::Float64=1e-5)::Int where T<:Number
    # check whether point lies on the boundary of poly
    for facet in poly.facets
        polygon = push!(map(v -> poly.verts[v], facet), poly.verts[facet[1]])
        if  inpolygon3d(point, polygon, tol) != 0
            return -1
        end
    end

    while true
        r = Ray(point, normalize!(randn(Float64, 3)))
        numIntersections = 0 # number of intersections of r and the facets of poly

        for facet in poly.facets
            polygon = push!(map(v -> poly.verts[v], facet), poly.verts[facet[1]])
            E = plane(polygon)

            try
                intersect(r, E)
            catch error
                continue
            end

            p = intersect(r, E)

            if inpolygon3d(p, polygon, tol) == -1
                break
            elseif inpolygon3d(p, polygon, tol) == 1
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



