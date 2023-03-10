using Plots

include("Polyhedron.jl")

# plot Polyhedron
function plotPolyhedron(poly::Polyhedron; fill_color = :red, kwargs...)
    p = plot()
    verts = poly.verts
    for facet in poly.facets
        xvals = append!(map(ind -> verts[ind][1], facet), verts[facet[1]][1])
        yvals = append!(map(ind -> verts[ind][2], facet), verts[facet[1]][2])
        zvals = append!(map(ind -> verts[ind][3], facet), verts[facet[1]][3])

        plot!(p, xvals, yvals, zvals, kwargs...)
    end
    return p
end