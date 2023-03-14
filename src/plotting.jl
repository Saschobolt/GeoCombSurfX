using PlotlyJS
using Colors

include("Polyhedron.jl")
include("decomposition.jl")

"""
Aux function to plot a polyhedron. Returns array of traces that can be handled by PlotlyJS.
"""
function tracePolyhedron(poly::Polyhedron; color::Color = RGB(0,0.9,1), text::Bool = false)
    polyTriang = triangulatePolyhedron(poly)

    facecolor = repeat([color], length(polyTriang.facets))

    # plot faces
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

    # plot edges
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

    # plot vertices
    mode = text ? "markers+text" : "markers"
    trace = scatter3d(
    	x = [v[1] for v in poly.verts],
        y = [v[2] for v in poly.verts],
        z = [v[3] for v in poly.verts],
        mode = mode,
        text = [string(v) for v in 1:length(poly.verts)],
        textposition = "bottom center",
        marker = attr(color = color, size = 4)
    )

    push!(traces, trace)

    return traces
end


function plotAssembly(assembly::Vector{Polyhedron}; text::Bool = false)
    colors = distinguishable_colors(length(assembly), RGB(0,0,0))

    plot(reverse(union([tracePolyhedron(poly; color = colors[i], text = text) for (i, poly) in enumerate(assembly)]...)), Layout(showlegend = false))
end


function plotPolyhedron(poly::Polyhedron; color::Color = RGB(0,0.9,1), text::Bool = false)
    plot(tracePolyhedron(poly; color = color, text = text), Layout(showlegend = false))
end