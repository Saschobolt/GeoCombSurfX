using PlotlyJS

include("../src/examples.jl")
include("../src/plotting.jl")
include("../src/decomposition.jl")

function create_plots(width = 600, height = 600, color_frame = RGB(1,0,0), color_blocks = RGB(0.5,0.5,0.5), color_links = RGB(0,0.25,1), showbackground = false, opacity = 1)
    # 4 prism
    p = plot(nprism(4), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "4prism.png", width = width, height = height)

    # 6 prism
    p = plot(nprism(6), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "6prism.png", width = width, height = height)

    # 8 prism
    p = plot(nprism(8), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "8prism.png", width = width, height = height)

    # merging 6 prism and 8 prism
    p = plot(merge(nprism(6), nprism(8), [[1,2,7,8]], [[1,2,9,10]]), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "merge1.png", width = width, height = height)

    p = plot(merge(nprism(6), nprism(8), [[1,2,7,8]], [[2,10,1,9]]), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "merge2.png", width = width, height = height)

    @info "Small blocks figures saved."

    # 3 hat
    hat = nprism(3)
    merge!(hat, nprism(3), [[3,1,6,4]], [[1,2,4,5]])
    merge!(hat, nprism(3), [[2,3,5,6]], [[1,2,4,5]])
    p = plot(flattenfacets(hat), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "hat.png", width = width, height = height)

    # tiblocks with 3 hat
    p = plot(flattenfacets(tiblock(6,3,2)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block632.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(8,3,2)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block832.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(6,3,3)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block633.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(8,3,4)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block834.png", width = width, height = height)

    @info "Blocks with 3hat figures saved."

    # assemblies with blocks above
    ass, frame = assembly1(6,3,2,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_6327.png", width = width, height = height)

    ass, frame = assembly1(8,3,2,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_8327.png", width = width, height = height)

    ass, frame = assembly2(6,3,2,7) 
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_6327.png", width = width, height = height)

    ass, frame = assembly2(8,3,2,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_8327.png", width = width, height = height)

    ass, frame, link = assembly3(6,4,3,1)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    colors[link] .= color_links
    p = plot(ass, colors = colors, width = width, height = height, showbackground = showbackground, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "assembly3_6431.png", width = width, height = height)

    ass, frame, link = assembly3(6,4,1,1)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    colors[link] .= color_links
    p = plot(ass, colors = colors, width = width, height = height, showbackground = showbackground, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "assembly3_6411.png", width = width, height = height)

    ass, frame = assembly1(6,3,3,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_6337.png", width = width, height = height)

    ass, frame = assembly1(8,3,4,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_8347.png", width = width, height = height)

    @info "Assemblies with 3hat figures saved."

    # tiblocks with 4 hat
    p = plot(flattenfacets(tiblock(6,4,2)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block642.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(8,4,2)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block842.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(6,4,3)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block643.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(8,4,4)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block844.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(9,4,3)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block943.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(12,4,4)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block1244.png", width = width, height = height)

    @info "Blocks with 4hat figures saved."

    # assemblies with blocks above
    ass, frame = assembly2(6,4,2,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_6427.png", width = width, height = height)

    ass, frame = assembly2(8,4,2,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_8427.png", width = width, height = height)

    ass, frame = assembly2(6,4,3,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_6437.png", width = width, height = height)

    ass, frame = assembly2(8,4,4,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_8447.png", width = width, height = height)

    ass, frame = assembly4(9,4,3,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly4_9437.png", width = width, height = height) 

    ass, frame = assembly4(12,4,4,7)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly4_12447.png", width = width, height = height) 

    @info "Assemblies with 4hat figures saved."

    # # Cube interlocking
    # ass, frame = cube_interlocking(10)
    # 
    # colors = [color_blocks for block in ass]
    # colors[frame] .= color_frame
    # p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    # PlotlyJS.savefig(p, "cube_interlocking.png", width = width, height = height) 

    # # Octahedra interlocking
    # ass, frame = octahedra_interlocking(10)
    # 
    # colors = [color_blocks for block in ass]
    # colors[frame] .= color_frame
    # p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    # PlotlyJS.savefig(p, "octahedra_interlocking.png", width = width, height = height) 

    # # Cube interlocking
    # ass, frame = tetrahedra_interlocking(10)
    # 
    # colors = [color_blocks for block in ass]
    # colors[frame] .= color_frame
    # p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    # PlotlyJS.savefig(p, "tetrahedra_interlocking.png", width = width, height = height)

    # edge-face contact
    block1 = deepcopy(Cube)
    block2 = deepcopy(Cube)
    aff = rigidmap(hcat([0,0,0], [1,0,0], [0,1,0], [0,0,1]), hcat([0.5,1,0], [0.5+1/sqrt(2), 1+1/sqrt(2), 0], [0.5-1/sqrt(2), 1+1/sqrt(2), 0], [0.5,1,1]))
    set_verts!(block2, aff.(get_verts(block2)))
    ass = [block1, block2]

    colors = [color_blocks for block in ass]
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = 0.55)
    for i in 1:2, j in 1:2
        add_trace!(
            p,
            scatter(
                x = [0.5, 0.5],
                y = [1, 1],
                z = [0, 1],
                line = attr(color = "red", width = 14),
                mode = "lines",
                type = "scatter3d"
            ),
            row = i, col = j
        )
    end

    PlotlyJS.savefig(p, "contact_edge_facet.png", width = width, height = height)

    # convex edge-edge contact
    block1 = deepcopy(Cube)
    block2 = deepcopy(Cube)
    aff = rigidmap(hcat([0,0,0], [1,0,0], [0,1,0], [0,0,1]), hcat([1,1,0], [2,1,0], [1,2,0], [1,1,1]))
    set_verts!(block2, aff.(get_verts(block2)))
    ass = [block1, block2]

    colors = [color_blocks for block in ass]
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = 0.55)
    for i in 1:2, j in 1:2
        add_trace!(
            p,
            scatter(
                x = [1, 1],
                y = [1, 1],
                z = [0, 1],
                line = attr(color = "red", width = 14),
                mode = "lines",
                type = "scatter3d"
            ),
            row = i, col = j
        )
    end

    PlotlyJS.savefig(p, "contact_convex_convex.png", width = width, height = height)

    # mixed edge-edge contact
    ass, frame = assembly1(6,3,2,2)
    colors = [color_blocks for block in ass]
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = 0.55)

    edges = [[10,11], [7,8], [8,9], [11,12]]
    block = ass[1]
    for i in 1:2, j in 1:2
        for e in edges
            add_trace!(
                p,
                scatter(
                    x = [get_verts(block)[e[1]][1], get_verts(block)[e[2]][1]],
                    y = [get_verts(block)[e[1]][2], get_verts(block)[e[2]][2]],
                    z = [get_verts(block)[e[1]][3], get_verts(block)[e[2]][3]],
                    line = attr(color = "red", width = 14),
                    mode = "lines",
                    type = "scatter3d"
                ),
                row = i, col = j
            )
        end
    end

    PlotlyJS.savefig(p, "contact_convex_concave.png", width = width, height = height)

    @info "Different contact points figures saved."
    @info "Finished."
end