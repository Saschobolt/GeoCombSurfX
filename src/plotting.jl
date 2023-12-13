using PlotlyJS
using Colors
import PlotlyJS.plot

include("Polyhedron.jl")
include("decomposition.jl")

function trace(f::AbstractEmbeddedGraph; 
                vertexcolors::AbstractVector{<:Color} = [RGB(0,0,0)], edgecolors::AbstractVector{<:Color} = [RGB(0,0,0)], labels::Bool = false,
                drawverts::Bool = false)
    if length(vertexcolors) == 1
        vertexcolors = repeat(vertexcolors, size(get_verts(f))[2])
    elseif length(vertexcolors) != size(get_verts(f))[2]
        error("Every vertex needs a color assigned to it.")
    end

    if length(edgecolors) == 1
        edgecolors = repeat(edgecolors, length(get_edges(f)))
    elseif length(edgecolors) != length(get_edges(f))
        error("Every edge needs a color assigned to it.")
    end
    traces = GenericTrace[]

    # plot edges
    for (i, e) in enumerate(get_edges(f))
        trace = scatter(
            x = get_verts(f)[1, e],
            y = get_verts(f)[2, e],
            z = get_verts(f)[3, e],

            line = attr(color = edgecolors[i]),
            mode = "lines",
            type = "scatter3d"
        )

        push!(traces, trace)
    end

    # plot vertices
    if drawverts
        mode = labels ? "markers+text" : "markers"
        trace = scatter3d(
            x = get_verts(f)[1, :],
            y = get_verts(f)[2, :],
            z = get_verts(f)[3, :],

            mode = mode,
            text = string.(collect(1:size(get_verts(f))[2])),
            textposition = "bottom center",
            marker = attr(color = vertexcolors)
        )

        push!(traces, trace)

        # for v in 1:size(get_verts(f))[2]
        #     trace = scatter3d(
        #         x = get_verts(f)[1, v],
        #         y = get_verts(f)[2, v],
        #         z = get_verts(f)[3, v],

        #         mode = mode,
        #         text = [string(v)], 
        #         textposition = "bottom center",
        #         marker = attr(color = vertexcolors[v])
        #     )

        #     push!(traces, trace)
        # end
    end
    
    return traces
end

function trace(poly::AbstractPolyhedron; 
            facetcolors::AbstractVector{<:Color} = [RGB(0,0.9,1)], opacity::Real = 0.5, kwargs...)
    polytriang = triangulate(poly)
    if length(facetcolors) == 1
        facetcolors = repeat(facetcolors, length(get_facets(polytriang)))
    elseif length(facetcolors) != length(get_facets(poly))
        error("Every facet needs a color assigned to it.")
    end

    mesh = mesh3d(
        x = get_verts(poly)[1,:],
        y = get_verts(poly)[2,:],
        z = get_verts(poly)[3,:],
        i = [triang[1] for triang in get_facets(polytriang)].-1,
        j = [triang[2] for triang in get_facets(polytriang)].-1,
        k = [triang[3] for triang in get_facets(polytriang)].-1,

        facecolor = facetcolors,
        opacity = opacity
    )

    traces = [mesh]
    append!(traces, trace(Framework(get_verts(poly), get_edges(poly)); kwargs...))

    return traces
end

function trace(assembly::AbstractVector{<:AbstractEmbeddedGraph}; 
        vertexcolors::AbstractVector{<:AbstractVector{<:Color}} = [[RGB(0,0,0)]], edgecolors::AbstractVector{<:AbstractVector{<:Color}} = [[RGB(0,0,0)]],
        facetcolors::AbstractVector{<:AbstractVector{<:Color}} = [[RGB(0,0.9,1)]],
        kwargs...)
    
    if length(vertexcolors) == 1
        vertexcolors = repeat(vertexcolors, length(assembly))
    elseif length(vertexcolors) != length(assembly)
        error("Every element of the assembly needs vertexcolors assigned to it.")
    end

    if length(edgecolors) == 1
        edgecolors = repeat(edgecolors, length(assembly))
    elseif length(edgecolors) != length(assembly)
        error("Every element of the assembly needs edgecolors assigned to it.")
    end

    if length(facetcolors) == 1
        facetcolors = repeat(facetcolors, length(assembly))
    elseif length(facetcolors) != length(assembly)
        error("Every element of the assembly needs facetcolors assigned to it.")
    end

    traces = vcat([trace(
            f; facetcolors = facetcolors[i], vertexcolors = vertexcolors[i], edgecolors = edgecolors[i], kwargs...
        ) for (i,f) in enumerate(assembly)]...)

    return traces
end

function plot(assembly::AbstractVector{<:AbstractEmbeddedGraph}; 
        subplots::Tuple{<:Integer, <:Integer} = (1,1), viewpoints::AbstractVecOrMat{<:Vector{<:Real}} = [[2.5,2.5,2.5]], width::Integer = 800, height::Integer = 800, showbackground::Bool = false, kwargs...)
    n_rows, n_cols = subplots
    n_subplots = n_rows * n_cols
    traces = trace(assembly; kwargs...)
    
    center_of_mass_verts = center_of_mass(hcat(get_verts.(assembly)...))
    p = make_subplots(rows=n_rows, cols=n_cols, specs=fill(Spec(kind="scene"), n_rows, n_cols), vertical_spacing = 0.05, horizontal_spacing = 0.05) # TODO: subplot matrix is only one row?
    # display(p)

    function subplot_layout(i)
        # set subplot camera eye and center
        subplot_eye = viewpoints[i]
        subplot_center = center_of_mass_verts

        # create scene
        scene = attr(
                xaxis_title = (showbackground ? "x" : ""),
                yaxis_title = (showbackground ? "y" : ""),
                zaxis_title = (showbackground ? "z" : ""),
                autosize = true, 
                aspectmode = "data",
                camera = attr(
                    eye = attr(x = subplot_eye[1], y = subplot_eye[2], z = subplot_eye[3])
                    # center = attr(x = subplot_center[1], y = subplot_center[2], z = subplot_center[3])
                ),
                xaxis = attr(
                    showbackground = showbackground,
                    showaxeslabels = showbackground,
                    showticklabels = showbackground
                ),
                yaxis = attr(
                    showbackground = showbackground,
                    showaxeslabels = showbackground,
                    showticklabels = showbackground
                ),
                zaxis = attr(
                    showbackground = showbackground,
                    showaxeslabels = showbackground,
                    showticklabels = showbackground
                )
            )   
        return scene
    end

    for i in 1:n_subplots
        # update the subplot layouts
        if i == 1
            relayout!(p, Dict(Symbol("scene") => subplot_layout(i)))
        else
            relayout!(p, Dict(Symbol("scene$(i)") => subplot_layout(i)))
        end

        # Compute subplot index in the matrix
        row_idx = div(i - 1, n_cols) + 1
        col_idx = (i - 1) % n_cols + 1

        # set subplot traces
        for trace in traces
            add_trace!(p, trace, row = row_idx, col = col_idx)
        end
    end
    
    relayout!(p, width = width, height = height, showlegend = false)
    return p
end

function plot(g::AbstractEmbeddedGraph; kwargs...)

    return plot([g]; kwargs...)
end
# """
# Aux function to plot a polyhedron. Returns array of traces that can be handled by PlotlyJS.
# """
# function trace_polyhedron(poly::AbstractPolyhedron; 
#             color::Color = RGB(0,0.9,1), labels::Bool = false, opacity::Real = 0.6, drawverts::Bool = true)
#     polyTriang = triangulate(poly)

#     facecolor = repeat([color], length(polyTriang.facets))

#     # plot faces
#     mesh = mesh3d(
#         x = [vert[1] for vert in polyTriang.verts],
#         y = [vert[2] for vert in polyTriang.verts],
#         z = [vert[3] for vert in polyTriang.verts],
#         i = [triang[1] for triang in polyTriang.facets].-1,
#         j = [triang[2] for triang in polyTriang.facets].-1,
#         k = [triang[3] for triang in polyTriang.facets].-1,
#         facecolor = facecolor,
#         opacity = opacity
#     )

#     traces = [mesh]

#     # plot edges
#     for edge in poly.edges
#         trace = scatter(
#             x = [poly.verts[v][1] for v in edge],
#             y = [poly.verts[v][2] for v in edge],
#             z = [poly.verts[v][3] for v in edge],
#             line = attr(color = "black"),
#             mode = "lines",
#             type = "scatter3d"
#         )

#         push!(traces, trace)
#     end

#     # plot vertices
#     if drawverts
#         mode = labels ? "markers+text" : "markers"
#         trace = scatter3d(
#             x = [v[1] for v in poly.verts],
#             y = [v[2] for v in poly.verts],
#             z = [v[3] for v in poly.verts],
#             mode = mode,
#             text = [string(v) for v in 1:length(poly.verts)],
#             textposition = "bottom center",
#             marker = attr(color = color, size = 4)
#         )
#     end
#     push!(traces, trace)

#     return traces
# end

# function plot(assembly::Vector{<:AbstractPolyhedron}; 
#               subfigures::Tuple{<:Integer, <:Integer} = (2,2), colors::Vector{<:Color} = [RGB(0,0.9,1)], labels::Bool = false, opacity::Real = 1, 
#               showbackground::Bool = true, drawverts::Bool = true, zoomfactor::Real = 1, viewpoint::Vector{<:Real} = [1,1,1], width::Int = 600, height::Int = 600)
#     if length(colors) == 1
#         colorvec = distinguishable_colors(length(assembly), colors[1])
#     else colorvec = colors
#     end

#     @assert length(colorvec) == length(assembly) "Only $(length(colors)) colors supplied for $(length(assembly)) blocks."

#     traces = reverse(union([trace_polyhedron(poly; color = colorvec[i], labels = labels, opacity = opacity, drawverts = drawverts) for (i, poly) in enumerate(assembly)]...))

#     scene(eyex = 1.25, eyey = 1.25, eyez = 1.25) = attr(
#         xaxis_title = (showbackground ? "x" : ""),
#         yaxis_title = (showbackground ? "y" : ""),
#         zaxis_title = (showbackground ? "z" : ""),
#         autosize = true, 
#         aspectmode = "data",
#         camera = attr(
#             eye = attr(x = eyex, y = eyey, z = eyez)
#         ),
#         xaxis = attr(
#             showbackground = showbackground,
#             showaxeslabels = showbackground,
#             showticklabels = false
#         ),
#         yaxis = attr(
#             showbackground = showbackground,
#             showaxeslabels = showbackground,
#             showticklabels = false
#         ),
#         zaxis = attr(
#             showbackground = showbackground,
#             showaxeslabels = showbackground,
#             showticklabels = false
#         )
#     )

#     fig = make_subplots(rows=subfigures[1], cols=subfigures[2], specs=fill(Spec(kind="scene"), 2, 2), vertical_spacing = 0.05, horizontal_spacing = 0.05)

#     if subfigures == (2,2)
#         relayout!(
#             fig, 
#             scene = scene((zoomfactor * 2.5 * normalize([1,1,1]))...),
#             scene2 = scene((zoomfactor * 2.5 * normalize([0,1,0.5]))...), # TODO: viewpoint als Argument einbauen.
#             scene3 = scene((zoomfactor * 2.5 * normalize([1,1,-1]))...),
#             scene4 = scene((zoomfactor * 2.5 * normalize([0,0,1]))...),
#             width = width,
#             height = height,
#             showlegend = false
#         )
#     elseif subfigures == (1,1)
#         relayout!(
#             fig,
#             scene = scene((zoomfactor * 2.5 * normalize(viewpoint))...),
#             showlegend = false
#         )
#     else error("Not implemented yet.") # TODO: Dynamische subfigures implementieren. Wie können scene-Attribute dynamisch gesetzt werden?
#     end

#     for i in 1:subfigures[1], j in 1:subfigures[2]
#         for trace in traces
#             add_trace!(
#                 fig,
#                 trace,
#                 row = i, col = j
#             )
#         end
#     end

#     fig
# end


# function plot(poly::AbstractPolyhedron; kwargs...)
#     plot([poly]; kwargs...)
# end