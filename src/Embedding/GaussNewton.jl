include("../Framework.jl")
include("../affine_geometry.jl")

"""
    gn_overconstrained(r::Function, J::Function, x0::AbstractVector{<:Real}; maxiter::Integer = 10000, tol::Real = 1e-12)

Gauss Newton method for overconstrained systems. 
"""
function gn_overconstrained(r::Function, J::Function, x0::AbstractVector{<:Real}; maxiter::Integer = 100000, tol::Real = 1e-12)
    x = copy(x0)

    for i in 1:maxiter
        div = J(x)
        del = - pinv(transpose(div) * div) * (transpose(div) * r(x))
        x += del
        if sqrt(sum(abs2.(del))) < tol
            break
        end
    end
    return x
end

"""
    gn_underonstrained(r::Function, J::Function, x0::AbstractVector{<:Real}; maxiter::Integer = 10000, tol::Real = 1e-12)

Gauss Newton method for underconstrained systems.
"""
function gn_underonstrained(r::Function, J::Function, x0::AbstractVector{<:Real}; maxiter::Integer = 100000, tol::Real = 1e-12)
    x = copy(x0)

    for _ in 1:maxiter
        div = J(x)
        del = - transpose(div) * pinv(div * transpose(div)) * r(x)
        x += del
        if sqrt(sum(abs2.(del))) < tol
            break
        end
    end
    return x
end

"""
    embed_gn(g::AbstractSimpleGraph, edgelengths::AbstractVector{<:Real}, init::AbstractVecOrMat(<:Real), d::Integer = 3; kwargs...)

Embed a graph into a d-dimensional space using Gauss-Newton method.
"""
function embed_gn(g::AbstractSimpleGraph, edgelengths::AbstractVector{<:Real}, init::AbstractVecOrMat{<:Real}, d::Integer = 3; kwargs...)
    if typeof(g) <: AbstractEmbeddedGraph # g is an embedded graph, thus get_verts returns the coordinate matrix
        n = size(get_verts(g))[2]
        n_vars = length(get_verts(g))
    else    # g is a combinatorial graph, thus get_verts returns an array [1, 2, ..., n] of vertex labels.
        n = length(get_verts(g))
    end
    n_vars = n*d
    n_constr = length(get_edges(g)) + Int(d*(d+1)/2 - max(0, d - length(get_verts(g))) * (d - length(get_verts(g)) + 1) / 2)

    # residual function
    function r(x, edgelengths)
        sol = zeros(n_constr)
        # distance constraints
        for (i, e) in enumerate(get_edges(g))
            sol[i] = sqdist(x[(e[1]-1)*d+1:e[1]*d], x[(e[2]-1)*d+1:e[2]*d]) - edgelengths[i]^2
        end

        next = length(get_edges(g)) + 1
        
        # constraints, that first vertex is [0,0,...], second is [x, 0,0,...], third is [x,y,0,0,...] and so forth
        for v in 1:min(length(get_verts(g)), d)
            for i in v:d
                sol[next] = x[(v-1)*d+i]^2
                next += 1
            end
        end
        return sol
    end

    # Jacobi Matrix of r at x
    function J(x, edgelengths)
        sol_tr = zeros(n_constr, n_vars)
        for (i,e) in enumerate(get_edges(g))
            sol_tr[(e[1]-1)*d+1:e[1]*d, i] = 2 * (x[(e[1]-1)*d+1:e[1]*d] - x[(e[2]-1)*d+1:e[2]*d])
            sol_tr[(e[2]-1)*d+1:e[2]*d, i] = -2 * (x[(e[1]-1)*d+1:e[1]*d] - x[(e[2]-1)*d+1:e[2]*d])
        end

        next = length(get_edges(g)) + 1

        for v in 1:min(length(get_verts(g)), d)
            for i in v:d
                sol_tr[(v-1)*d+i, next] = 2 * x[(v-1)*d+i]
                next += 1
            end
        end
        return transpose(sol_tr)
    end

    if n_constr >= n_vars
        coords = reshape(gn_overconstrained(x -> r(x, edgelengths), x -> J(x, edge_lengths), reshape(init, n_vars); kwargs...), d, :)
    elseif n_constr < n_vars
        coords = reshape(gn_underonstrained(x -> r(x, edgelengths), x -> J(x, edge_lengths), reshape(init, n_vars); kwargs...), d, :)
    end

    return Framework(coords, get_edges(g))
end