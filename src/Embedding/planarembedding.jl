include("../SimplicialSurface.jl")

using JuMP
using HiGHS

function planar_embedding_angles(surf::AbstractColoredSimplicialSurface, boundangles_dict::Dict{<:Integer,<:Real}=Dict{Int,Float64}())
    cols = colors(surf)
    outer = unique(vcat(boundary(surf)...))
    inner = setdiff(surf.verts, outer)
    @assert issubset(keys(boundangles_dict), outer) "The keys of boundangles_dict need to be boundary vertices of surf, but got $(collect(keys(boundangles_dict))) instead."

    model = Model(HiGHS.Optimizer)

    # Variables
    @variable(model, 0 <= anglesum_at_vertex[surf.verts] <= 2 * π)
    @variable(model, 0 <= angle_between_cols[cols, cols] <= 2 * π)

    # Constraints
    # orientation doesn't matter for angle sum
    for i in eachindex(cols)
        for j in i+1:length(cols)
            @constraint(model, angle_between_cols[cols[i], cols[j]] == angle_between_cols[cols[j], cols[i]])
        end
    end

    # angle sum at inner vertices is 2π
    for v in inner
        @constraint(model, anglesum_at_vertex[v] == 2 * π)
    end

    # angle sum at boundary vertices is given
    for (v, anglesum) in boundangles_dict
        @constraint(model, anglesum_at_vertex[v] == anglesum)
    end

    # angle sum at a vertex is the sum of the angles between the colored legs of the incident faces
    for v in surf.verts
        faces = incfacets(surf, v)
        legs = [[color(surf, [v, setdiff(f, v)[1]]), color(surf, [v, setdiff(f, v)[2]])] for f in faces]
        @constraint(model, sum(angle_between_cols[leg...] for leg in legs) == anglesum_at_vertex[v])
    end

    # inner angle sum of triangles is π
    for f in surf.facets
        legs = [[color(surf, [v, setdiff(f, v)[1]]), color(surf, [v, setdiff(f, v)[2]])] for v in f]
        @constraint(model, sum(angle_between_cols[leg...] for leg in legs) == π)
    end

    # Objective: Don't want very thin triangles. So maximize the minimum angle between colored legs by introducing a new variable smaller than all angles and maximizing it.
    @variable(model, 0 <= minangle <= 2 * π)
    for i in eachindex(cols)
        for j in i+1:length(cols)
            @constraint(model, minangle <= angle_between_cols[cols[i], cols[j]])
        end
    end
    @objective(model, Max, minangle)

    optimize!(model)

    if !is_solved_and_feasible(model)
        @info "There exist no feasible triangulations as corresponding LP doesn't have a feasible solution. Check solution_summary to recieve more details."
    else
        @info "There exist a feasible planar triangulation with similarity classes defined by the edge coloring. The angles between colors are:\n
        $(value.(angle_between_cols))\n
        The angle sums at the vertices are \n
        $(value.(anglesum_at_vertex))."
    end

    return model
end