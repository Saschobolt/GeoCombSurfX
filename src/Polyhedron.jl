using GeometryBasics

struct Polyhedron
    verts::Array{1, Array{1, Float32}} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Array{1, Tuple{Int}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Array{1, Array{1, Int}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
end

# plot Polyhedron

# triangulate surface of Polyhedron

# combinatorics