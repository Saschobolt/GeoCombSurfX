using Plots
using PolygonOps
using LinearAlgebra


struct Polyhedron
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


# plot Polyhedron
function plotPolyhedron!(p::Plots.Plot, poly::Polyhedron; fill_color = :red, kwargs...)
    for facet in poly.facets
        xvals = append!(map(ind -> verts[ind][1], facet), verts[facet[1]][1])
        yvals = append!(map(ind -> verts[ind][2], facet), verts[facet[1]][2])
        zvals = append!(map(ind -> verts[ind][3], facet), verts[facet[1]][3])

        plot!(p, xvals, yvals, zvals, kwargs...)
    end
    return p
end

"""
A::Matrix
cond::Float64
Returns an orthonormal basis for the columnspace of the matrix A using svd. Singular values with abs < cond are treated as 0.
"""
function colspace(A::Matrix{T}, cond::Float64 = 1e-8)::Vector{Vector{Float64}} where T<:Number
    F = svd(A)
    return [c[:] for c in eachcol(F.U[:, findall(>(cond), abs.(F.S))])]
end


"""
A::Matrix
cond::Float64
Returns an orthonormal basis for the nullspace of the matrix A using svd. Singular values with abs < cond are treated as 0.
"""
function nullspace(A::Matrix{T}, cond::Float64 = 1e-8)::Vector{Vector{Float64}} where T<:Number
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

    A = hcat(union(plane.vectors, ray.vector)...)
    if length(colspace(A)) < 3
        throw(ErrorException("ray and plane are parallel."))
    end
    
    coeffs = A\(-vp-vr)
    if coeffs[3] < 0
        throw(ErrorException("ray and plane don't intersect"))
    end

    return vr + coeffs[3] * ray.vector
end


"""
returns 1 if the 3d point p lies in the polygon poly embedded into R^3, -1 if it lies on its boundary and 0 otherwise
"""
function inpolygon3d(p::AbstractVector, poly::AbstractVector, tol::Float64=1e-8)
    E = plane(poly)
    point = p

    # display(hcat(E.vectors...))
    # print((hcat(E.vectors...) * (hcat(E.vectors...)\(point-E.point))) - (point-E.point))

    if max(abs.((hcat(E.vectors...) * (hcat(E.vectors...)\(point-E.point))) - (point-E.point))...) < tol
        return inpolygon(point, poly) # inpolygon projects onto the xy plane and decides the problem. If p and poly lie in the same plane, then this is enough to decide the problem
    else
        return 0 # p and poly don't lie in the same plane
    end
end

# triangulate surface of Polyhedron
"""
returns a polyhedron containing the vertices and edges of poly such that every facet is triangular.
"""
function triangulatePolyhedron(poly::Polyhedron)::Polyhedron
    newVerts = poly.verts
    newEdges = poly.edges
    newFacets = []


    for facet in poly.facets
        subfacet = facet

        while length(subfacet) > 3
            if inpolygon((poly.verts[subfacet[length(subfacet)]] + poly.verts[subfacet[2]]) / 2, append!(map(i -> poly.verts[i], subfacet), [poly.verts[subfacet[1]]])) == 1
                append!(newEdges, [(subfacet[length(subfacet)], subfacet[2])])
                append!(newFacets, [[subfacet[length(subfacet)], subfacet[1], subfacet[2]]])
                subfacet = subfacet[2:length(subfacet)]
                continue
            elseif inpolygon((poly.verts[subfacet[length(subfacet) - 1]] + poly.verts[subfacet[1]]) / 2, append!(map(i -> poly.verts[i], subfacet), [poly.verts[subfacet[1]]])) == 1
                append!(newEdges, [(subfacet[length(subfacet) - 1], subfacet[1])])
                append!(newFacets, [[subfacet[length(subfacet) - 1], subfacet[length(subfacet)], subfacet[1]]])
                subfacet = subfacet[1:length(subfacet) - 1]
                continue
            else
                for ind in 2:length(subfacet) - 1
                    if inpolygon((poly.verts[subfacet[ind - 1]] + poly.verts[subfacet[ind + 1]]) / 2, append!(map(i -> poly.verts[i], subfacet), [poly.verts[subfacet[1]]])) == 1
                        append!(newEdges, [(subfacet[ind - 1], subfacet[ind + 1])])
                        append!(newFacets, [[subfacet[ind - 1], subfacet[ind], subfacet[ind + 1]]])
                        subfacet = subfacet[1:end .!= ind]
                        break
                    end
                end
            end
        end

        append!(newFacets, [subfacet])
    end

    return Polyhedron(newVerts, newEdges, newFacets)
end


# decide whether point is inside of polyhedron


"""
    Randomized algorithm to check whether a point is contained in a polyhedron.
"""
function inpolyhedron(point::Vector{T}, poly::Polyhedron)::Int where T<:Number
    # check whether point lies on the boundary of poly
    for facet in poly.facets
        polygon = append!(map(v -> poly.verts[v], facet), [poly.verts[facet[1]]])
        if  inpolygon3d(point, polygon) != 0
            return -1
        end
    end

    while true
        r = Ray(point, normalize!(randn(Float64, 3)))
        numIntersections = 0 # number of intersections of r and the facets of poly

        for facet in poly.facets
            polygon = append!(map(v -> poly.verts[v], facet), [poly.verts[facet[1]]])
            E = plane(polygon)
            try
                p = intersect(r, E)
            catch error
                continue
            end

            if inpolygon3d(p, polygon) == -1
                break
            elseif inpolygon3d(p, polygon) == 1
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
