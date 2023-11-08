using PlotlyJS
using Colors
import PlotlyJS.plot

include("Polyhedron.jl")
include("decomposition.jl")

"""
Aux function to plot a polyhedron. Returns array of traces that can be handled by PlotlyJS.
"""
function trace_polyhedron(poly::AbstractPolyhedron; 
            color::Color = RGB(0,0.9,1), labels::Bool = false, opacity::Real = 0.6, drawverts::Bool = true)
    polyTriang = triangulate(poly)

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
        opacity = opacity
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
    if drawverts
        mode = labels ? "markers+text" : "markers"
        trace = scatter3d(
            x = [v[1] for v in poly.verts],
            y = [v[2] for v in poly.verts],
            z = [v[3] for v in poly.verts],
            mode = mode,
            text = [string(v) for v in 1:length(poly.verts)],
            textposition = "bottom center",
            marker = attr(color = color, size = 4)
        )
    end
    push!(traces, trace)

    return traces
end

function plot(assembly::Vector{<:AbstractPolyhedron}; 
              subfigures::Tuple{<:Integer, <:Integer} = (2,2), colors::Vector{<:Color} = [RGB(0,0.9,1)], labels::Bool = false, opacity::Real = 1, 
              showbackground::Bool = true, drawverts::Bool = true, zoomfactor::Real = 1, viewpoint::Vector{<:Real} = [1,1,1], width::Int = 600, height::Int = 600)
    if length(colors) == 1
        colorvec = distinguishable_colors(length(assembly), colors[1])
    else colorvec = colors
    end

    @assert length(colorvec) == length(assembly) "Only $(length(colors)) colors supplied for $(length(assembly)) blocks."

    traces = reverse(union([trace_polyhedron(poly; color = colorvec[i], labels = labels, opacity = opacity, drawverts = drawverts) for (i, poly) in enumerate(assembly)]...))

    scene(eyex = 1.25, eyey = 1.25, eyez = 1.25) = attr(
        xaxis_title = (showbackground ? "x" : ""),
        yaxis_title = (showbackground ? "y" : ""),
        zaxis_title = (showbackground ? "z" : ""),
        autosize = true, 
        aspectmode = "data",
        camera = attr(
            eye = attr(x = eyex, y = eyey, z = eyez)
        ),
        xaxis = attr(
            showbackground = showbackground,
            showaxeslabels = showbackground,
            showticklabels = false
        ),
        yaxis = attr(
            showbackground = showbackground,
            showaxeslabels = showbackground,
            showticklabels = false
        ),
        zaxis = attr(
            showbackground = showbackground,
            showaxeslabels = showbackground,
            showticklabels = false
        )
    )

    fig = make_subplots(rows=subfigures[1], cols=subfigures[2], specs=fill(Spec(kind="scene"), 2, 2), vertical_spacing = 0.05, horizontal_spacing = 0.05)

    if subfigures == (2,2)
        relayout!(
            fig, 
            scene = scene((zoomfactor * 2.5 * normalize([1,1,1]))...),
            scene2 = scene((zoomfactor * 2.5 * normalize([0,1,0.5]))...), # TODO: viewpoint als Argument einbauen.
            scene3 = scene((zoomfactor * 2.5 * normalize([1,1,-1]))...),
            scene4 = scene((zoomfactor * 2.5 * normalize([0,0,1]))...),
            width = width,
            height = height,
            showlegend = false
        )
    elseif subfigures == (1,1)
        relayout!(
            fig,
            scene = scene((zoomfactor * 2.5 * normalize(viewpoint))...),
            showlegend = false
        )
    else error("Not implemented yet.") # TODO: Dynamische subfigures implementieren. Wie können scene-Attribute dynamisch gesetzt werden?
    end

    for i in 1:subfigures[1], j in 1:subfigures[2]
        for trace in traces
            add_trace!(
                fig,
                trace,
                row = i, col = j
            )
        end
    end

    fig
end


function plot(poly::AbstractPolyhedron; kwargs...)
    plot([poly]; kwargs...)
end