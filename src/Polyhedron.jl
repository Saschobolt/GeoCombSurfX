using Plots

struct Polyhedron
    verts::Vector{Vector{Float32}} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{Tuple{Int, Int}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{Vector{Int}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
end

# plot Polyhedron
function plot_polyhedron!(p::Plots.Plot, poly::Polyhedron; fill_color = :red, kwargs...)
    for facet in poly.facets
        xvals = append!(map(ind -> verts[ind][1], facet), verts[facet[1]][1])
        yvals = append!(map(ind -> verts[ind][2], facet), verts[facet[1]][2])
        zvals = append!(map(ind -> verts[ind][3], facet), verts[facet[1]][3])

        plot!(p, xvals, yvals, zvals ,fill = (0, 0.5, fill_color), kwargs...)
    end
    return p
end


# triangulate surface of Polyhedron

# combinatorics