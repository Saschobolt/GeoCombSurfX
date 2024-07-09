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

# Domino partitions
import Base.display
# import Base.isequal
# import Base.==
# import Base.sum
mutable struct Piece
    n::Int # Integer representing the number of the piece
    in::Int # Integer representing the color of the left side
    out::Int # Integer representing the color of the right side
    twin::Union{Piece,Nothing}

    function Piece(n, in, out)
        self = new(n, in, out, nothing)
        twin = new(n, out, in, self)
        self.twin = twin
        return self
    end
end

function isequal(p1::Piece, p2::Piece)
    return p1.n == p2.n && p1.in == p2.in && p1.out == p2.out
end

function ==(p1::Piece, p2::Piece)
    return isequal(p1, p2)
end

function display(p::Piece)
    println("$(p.in)    $(p.out)\n $(p.n)")
end

mutable struct domino_partition
    pieces::Vector{Piece}

    function domino_partition(pieces::Vector{Piece})
        partition = new(pieces)
        for i in 1:length(pieces)-1
            if pieces[i].out != pieces[i+1].in
                error("In a domino partition, the out of a piece should be the in of the next piece, but out of the piece $(i) is $(pieces[i].out) and in of the piece $(i+1) is $(pieces[i+1].in).")
            end
        end

        return partition
    end
end

function sum(p::domino_partition)
    return Base.sum(piece.n for piece in p.pieces)
end

function in(p::domino_partition)
    return p.pieces[1].in
end

function out(p::domino_partition)
    return p.pieces[end].out
end

function isequal(p1::domino_partition, p2::domino_partition)
    return p1.pieces == p2.pieces
end

function ==(p1::domino_partition, p2::domino_partition)
    return isequal(p1, p2)
end

function domino_partition(n::Int, in::Int, out::Int, pieces::Vector{Piece})
    # calculate a domino partition of n consisting of the pieces in pieces, such that the first piece has in as its in and the last piece has out as its out. 
    
end