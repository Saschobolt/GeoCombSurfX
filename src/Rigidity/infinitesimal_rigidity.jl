include("../Framework.jl")
include("../affine_geometry.jl")

"""
    rigidity_matrix(f::AbstractEmbeddedGraph)

Construct the rigidity matrix of the embedded graph f.
"""
function rigidity_matrix(f::AbstractEmbeddedGraph)
    d, n = size(get_verts(f))
    r = zeros((length(get_edges(f)), d*n))
    coords = get_verts(f)
    for (i, e) in enumerate(get_edges(f))
        v = e[1]
        w = e[2]
        r[i, (v-1)*d+1:v*d] = coords[:, v] - coords[:, w]
        r[i, (w-1)*d+1:w*d] = coords[:, w] - coords[:, v]
    end

    return r
end

"""
    basis_inf_motions(f::AbstractEmbeddedGraph)

Calculate a basis for the space of infinitesimal motions of the embedded graph f.
"""
function basis_inf_motions(f::AbstractEmbeddedGraph; atol = 1e-10)
    return nullspace(rigidity_matrix(f), atol = atol)
end

function is_infrigid(f::AbstractEmbeddedGraph)
    d = dimension(f)
    return size(basis_inf_motions(f))[2] == binomial(d+1, 2)
end

"""
    basis_inf_flex(f::AbstractEmbeddedGraph)

Calculate a basis for a complement of the space of infinitesimal rigid motions of the embedded graph f.
"""
function basis_inf_flex(f::AbstractEmbeddedGraph)
    # TODO!
end

"""
    index(f::AbstractEmbeddedGraph)

Calculate the index of the embedded graph f.
"""
function index(f::AbstractEmbeddedGraph)
    d, n = size(get_verts(f))
    return length(get_edges(f)) - d*n + binomial(d+1, 2)
end