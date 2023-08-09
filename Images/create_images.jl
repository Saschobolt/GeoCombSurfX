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

    # assemblies with blocks above
    ass, frame = assembly1(6,3,2,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_6327.png", width = width, height = height)

    ass, frame = assembly1(8,3,2,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_8327.png", width = width, height = height)

    ass, frame = assembly2(6,3,2,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_6327.png", width = width, height = height)

    ass, frame = assembly2(8,3,2,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_8327.png", width = width, height = height)

    ass, frame, link = assembly3(6,4,3,1)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    colors[link] .= color_links
    p = plot(ass, colors = colors, width = width, height = height, showbackground = showbackground, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "assembly3_6431.png", width = width, height = height)

    ass, frame, link = assembly3(6,4,1,1)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    colors[link] .= color_links
    p = plot(ass, colors = colors, width = width, height = height, showbackground = showbackground, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "assembly3_6411.png", width = width, height = height)

    ass, frame = assembly1(6,3,3,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_6337.png", width = width, height = height)

    ass, frame = assembly1(8,3,4,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_8347.png", width = width, height = height)

    # tiblocks with 4 hat
    p = plot(flattenfacets(tiblock(6,4,2)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block642.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(8,4,2)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block842.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(6,4,3)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block643.png", width = width, height = height)

    p = plot(flattenfacets(tiblock(8,4,4)), showbackground = showbackground, width = width, height = height, drawverts = false, opacity = opacity)
    PlotlyJS.savefig(p, "block844.png", width = width, height = height)

    # assemblies with blocks above
    ass, frame = assembly2(6,4,2,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_6427.png", width = width, height = height)

    ass, frame = assembly2(8,4,2,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly2_8427.png", width = width, height = height)

    ass, frame = assembly1(6,4,3,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_6437.png", width = width, height = height)

    ass, frame = assembly1(8,4,4,7)
    ass = flattenfacets.(ass)
    colors = [color_blocks for block in ass]
    colors[frame] .= color_frame
    p = plot(ass, colors = colors, width = width, height = height, drawverts = false, showbackground = showbackground, opacity = opacity)
    PlotlyJS.savefig(p, "assembly1_8447.png", width = width, height = height)
end