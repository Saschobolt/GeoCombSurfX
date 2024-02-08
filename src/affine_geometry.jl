using LinearAlgebra

# """
# A::Matrix
# cond::Float64
# Returns an orthonormal basis for the nullspace of the matrix A using QR decomposition. Singular values with abs < cond are treated as 0.
# """
# function nullspace(A::Matrix{<:Real}, atol::Real = 1e-8)
#     Q, R = qr(A)
#     basis = Q[:, findall(x -> abs(x) < atol, diag(R))]
#     basis = [basis[:, i] for i in 1:size(basis)[2]]
#     return basis
# end

"""
    indcols_indices(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)


calculate the indices of the columns of A forming a maximal linear independent subset.
"""
function indcols_indices(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)
    Q, R, perm = qr(A, ColumnNorm())
    indices = perm[findall(val -> abs(val) > atol, diag(R))]
    return indices
end


"""
    indcols(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)

calculate maximal linear independent subset of the columns of A.
"""
function indcols(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)
    indices = indcols_indices(A, atol = atol)
    return A[:, indices]
end

"""
    colspace(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)


Returns an orthonormal basis for the columnspace of the matrix A using QR decomposition. Singular values with abs < atol are treated as 0.
"""
function colspace(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)
    Q, R = qr(indcols(A))
    return Q[:, 1:size(R)[1]]
end

"""
    affinebasis_indices(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)

Find the column indices of the Matrix A that form an affine basis of the affine space spanned by the columns of A.
Real values < atol are considered 0.
"""
function affinebasis_indices(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)
    M = vcat(A, transpose(repeat([1], size(A)[2])))
    return indcols_indices(M, atol = atol)
end

"""
    affinebasis_indices(s::Vector{<:Vector{<:Real}})

Find the indices of the entries in s that form an affine basis of the affine space spanned by the entries of s.
Real values < atol are considered 0.
"""
function affinebasis_indices(s::Vector{<:Vector{<:Real}}; atol::Real = 1e-8)
    return affinebasis_indices(hcat(s...), atol = atol)
end

"""
    affinebasis(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)

Find an affine basis of the affine space spanned by the points columns of s.
Real values < atol are considered zero.
"""
function affinebasis(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)
    indices = affinebasis_indices(A, atol = atol)
    return A[:,indices]
end

"""
    affinebasis(s::Vector{<:Vector{<:Real}})

Find an affine basis of the affine space spanned by the points in s.
Real values < atol are considered zero.
"""
function affinebasis(s::Vector{<:Vector{<:Real}}; atol::Real = 1e-8)
    return affinebasis(hcat(s...), atol = atol)
end

"""
    affinedim(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)

Affine dimension of the affine space spanned by the columns of A.
Real values < atol are considered zero.
"""
function affinedim(A::AbstractMatrix{<:Real}; atol::Real = 1e-8)
    return rank(vcat(A, transpose(repeat([1], size(A)[2]))), atol = atol) - 1
end

"""
    affinedim(s::Vector{<:Vector{<:Real}}; atol::Real = 1e-8)

Affine dimension of the affine space spanned by the entries of s.
Real values < atol are considered zero.
"""
function affinedim(s::Vector{<:Vector{<:Real}}; atol::Real = 1e-8)
    return affinedim(hcat(s...), atol = atol)
end

"""
    affinemap(preim::AbstractMatrix{<:Real}, im::AbstractMatrix{<:Real}; atol = 1e-8)

TBW
"""
function affinemap(preim::AbstractMatrix{<:Real}, im::AbstractMatrix{<:Real}; atol = 1e-8)
    @assert size(preim)[2] == size(im)[2] "Number of preim elements ($(size(preim)[2]) and image elements ($(size(im)[2])) need to match.)"

    d_pre = size(preim)[1] # dimension of underlying space of preimage
    @assert affinedim(preim) == d_pre "preim needs to contain an affine basis, but span of points has affine dimension $(affinedim(preim)) < $(d_pre)"
    d_im = size(im)[1] # dimension of underlying space of image

    basisind = affinebasis_indices(preim, atol = atol)
    
    A = preim[:, basisind]
    A = vcat(A, transpose(repeat([1], d_pre + 1))) # embed preimage into higher dimensional space
    b = im[:, basisind]
    b = vcat(b, transpose(repeat([1], d_pre + 1))) # embed image into higher dimensional space

    function aff(x::AbstractVecOrMat{<:Real})
        if typeof(x) <: AbstractVector
            y = vcat(x, [1])
        elseif typeof(x) <: AbstractMatrix
            y = vcat(x, ones(1, size(x)[2]))
        end

        M = b * inv(A)
        sol = (M * y)[1:(end-1), :]
        if size(sol)[2] == 1
            return sol[:, 1]
        end

        return (M * y)[1:(end-1), :]
    end

    return aff
end


"""
    affinemap(preim::Vector{<:Vector{<:Real}}, im::Vector{<:Vector{<:Real}})

TBW
"""
function affinemap(preim::Vector{<:Vector{<:Real}}, im::Vector{<:Vector{<:Real}}; atol = 1e-8)
    return affinemap(hcat(preim...), hcat(im...), atol = atol)
end

"""
    rigidmap(preim::Vector{<:Vector{<:Real}}, im::Vector{<:Vector{<:Real}})

TBW
"""
function rigidmap(preim::Matrix{<:Real}, im::Matrix{<:Real}; atol::Real = 1e-8)
    @assert size(preim)[2] == size(im)[2] "Number of preim elements ($(size(preim)[2]) and image elements ($(size(im)[2])) need to match.)"

    basisind = affinebasis_indices(preim)
    preimbasis = preim[:, basisind]
    imbasis = im[:, basisind]

    for i in 1:size(preimbasis)[2]
        for j in (i+1):size(preimbasis)[2]
            @assert abs(dist(preimbasis[:,i], preimbasis[:,j]) - dist(imbasis[:,i], imbasis[:,j])) < atol "Distance between preimage and image points need needs to be identical, but the distance between the points $(i) and $(j) is $(dist(preimbasis[:,i], preimbasis[:,j])) in the perimage and $(dist(imbasis[:,i], imbasis[:,j])) in the image."
        end
    end

    return affinemap(preim, im, atol=atol)
end

"""
    rigidmap(preim::Vector{<:Vector{<:Real}}, im::Vector{<:Vector{<:Real}})

TBW
"""
function rigidmap(preim::Vector{<:Vector{<:Real}}, im::Vector{<:Vector{<:Real}}; atol::Real = 1e-8)
    return rigidmap(hcat(preim...), hcat(im...), atol = atol)
end

struct Ray{T<:Real}
    point::Vector{T}
    vector::Vector{T}

    function Ray(point::AbstractVector{<:Real}, vector::AbstractVector{<:Real})
        return new{Float64}(point, vector)
    end
end

struct Plane{T<:Real}
    point::Vector{T} 
    vectors::Matrix{T}

    function Plane(point::AbstractVector{<:Real}, vectors::AbstractMatrix{<:Real})
        return new{Float64}(point, vectors)
    end
end

function Plane(points::AbstractMatrix{<:Real}; atol::Real=1e-8)
    @assert affinedim(points) == 2 "points doesn't span a plane."
    basis = affinebasis(points, atol = atol)
    
    v = basis[:, 1]
    A = basis[:, 2:end] - repeat(v, 1, size(basis)[2] - 1)
    return Plane(v, A)
end

"""
    Plane(points::Vector{<:Vector{<:Real}}; atol::Real=1e-8)

Calculate the affine plane in which points lie.
"""
function Plane(points::Vector{<:Vector{<:Real}}; atol::Real=1e-8)
    return Plane(hcat(points...))
end

"""
    normalvec(plane::Plane)

Calculate the normalized normal vector of a plane.
"""
function normalvec(plane::Plane)
    return normalize(cross(plane.vectors[:, 1], plane.vectors[:, 2]))
end

"""
    normalvec(polygon::Vector{<:Vector{<:Real}})

Calculate the normal vector of the plane that is spanned by points.
"""
function normalvec(points::Vector{<:Vector{<:Real}})
    return normalvec(Plane(points))
end

function normalvec(points::AbstractMatrix{<:Real})
    return normalvec([points[:,i] for i in 1:size(points)[2]])
end

"""
returns the intersection between a ray and a plane if it exists.
"""
function intersect(ray::Ray, plane::Plane; atol::Real = 1e-8)
    vr = ray.point
    vp = plane.point

    A = hcat(plane.vectors, -ray.vector)
    if rank(A) < 3
        throw(ErrorException("ray and plane are parallel."))
    end
    
    coeffs = A\(vr - vp)
    if coeffs[3] < 0
        throw(ErrorException("ray and plane don't intersect"))
    end

    return vr + coeffs[3] * ray.vector
end

function dist(v::Vector{<:Real}, w::Vector{<:Real})
    return norm(v-w)
end

function sqdist(v::Vector{<:Real}, w::Vector{<:Real})
    return sum((v-w).^2)
end

"""
    signedangle3d_right(v::Vector{<:Real}, w::Vector{<:Real}, n::Vector{<:Real}; atol = 1e-12)

Calculate the signed angle of the right handed rotation from the vector v to w with regard to the plane normal vector n. Real values < atol are considered 0. 
"""
function signedangle3d_right(v::Vector{<:Real}, w::Vector{<:Real}, n::Vector{<:Real}; atol = 1e-12)
    # https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
    @assert length(v) == 3 "Inputs need to be real vectors of length 3."
    @assert length(w) == 3 "Inputs need to be real vectors of length 3."
    @assert abs(dot(v,n)) < atol "v and n are not perpendicular. abs(dot(v,n)) = $(abs(dot(v,n)))"
    @assert abs(dot(w,n)) < atol "w and n are not perpendicular. abs(dot(w,n)) = $(abs(dot(w,n)))"
    return atan(dot(cross(v,w), n), dot(v,w))
end

"""
    signedangle3d_left(v::Vector{<:Real}, w::Vector{<:Real}, n::Vector{<:Real}; atol = 1e-12)

Calculate the signed angle of the left handed rotation from the vector v to w with regard to the plane normal vector n. Real values < atol are considered 0. 
"""
function signedangle3d_left(v::Vector{<:Real}, w::Vector{<:Real}, n::Vector{<:Real}; atol = 1e-12)
    # https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
    return -signedangle3d_right(v,w,n)
end

"""
    center_of_mass(points::Vector{<:Vector{<:Real}})

Calculate the center of a set of points.
"""
function center_of_mass(points::Vector{<:Vector{<:Real}})
    return sum(points) / length(points)
end

function center_of_mass(points::AbstractMatrix{<:Real})
    return center_of_mass([points[:,i] for i in 1:size(points)[2]])
end