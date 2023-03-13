using PlotlyJS
using Colors

include("Polyhedron.jl")
include("decomposition.jl")

"""
Aux function to plot a polyhedron. Returns array of traces that can be handled by PlotlyJS.
"""
function tracePolyhedron(poly::Polyhedron, color::Color = RGB(0,0.6,1))
    polyTriang = triangulatePolyhedron(poly)

    facecolor = repeat([color], length(polyTriang.facets))

    mesh = mesh3d(
        x = [vert[1] for vert in polyTriang.verts],
        y = [vert[2] for vert in polyTriang.verts],
        z = [vert[3] for vert in polyTriang.verts],
        i = [triang[1] for triang in polyTriang.facets].-1,
        j = [triang[2] for triang in polyTriang.facets].-1,
        k = [triang[3] for triang in polyTriang.facets].-1,
        facecolor = facecolor,
        opacity = 0.6
    )

    traces = [mesh]

    for edge in poly.edges
        trace = scatter(
            x = [poly.verts[v][1] for v in edge],
            y = [poly.verts[v][2] for v in edge],
            z = [poly.verts[v][3] for v in edge],
            line = attr(color = "black"),
            mode = "lines",
            type = "scatter3d"
        )

        push!(traces, trace)
    end

    return traces
end


function plotAssembly(assembly::Vector{Polyhedron})
    colors = distinguishable_colors(length(assembly))

    plot(union([tracePolyhedron(poly, colors[i]) for (i, poly) in enumerate(assembly)]...))
end


function plotPolyhedron(poly::Polyhedron; color::Color = RGB(0.2,0.2,1))
    plot(tracePolyhedron(poly, color))
end