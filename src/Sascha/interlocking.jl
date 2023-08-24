using JuMP
import Ipopt
using Polyhedra
import MultiObjectiveAlgorithms as MOA

include("../Polyhedron.jl")
include("../affine_geometry.jl")
include("../decomposition.jl")

"""
    contactfacets(poly1::Polyhedron, poly2::Polyhedron; atol = 1e-8)

TBW
"""
function contactfacets(poly1::Polyhedron, poly2::Polyhedron; atol = 1e-8)
    coords_facets1 = [get_verts(poly1)[facet] for facet in get_facets(poly1)]
    coords_facets2 = [get_verts(poly2)[facet] for facet in get_facets(poly2)]

    poly_facets1 = [Polyhedra.polyhedron(vrep(transpose(hcat(coords...)))) for coords in coords_facets1]
    poly_facets2 = [Polyhedra.polyhedron(vrep(transpose(hcat(coords...)))) for coords in coords_facets2]

    intersections = Matrix{Vector{Vector{<:Real}}}([[vrep(Polyhedra.intersect(face1, face2)).V[k,:] for k in 1:size(vrep(Polyhedra.intersect(face1, face2)).V)[1]] for face1 in poly_facets1, face2 in poly_facets2])

    return intersections
end


function contactfacets(assembly::Vector{<:Polyhedron}; atol = 1e-8)
    contacts = Dict()

    for i in 1:length(assembly)
        for j in i+1:length(assembly)
            intersections = contactfacets(assembly[i], assembly[j])
            if any(length.(intersections) .> 0)
                get!(contacts, [i,j], intersections)
            end
        end
    end 

    return contacts
end


"""
    titest(assembly::Vector{<:Polyhedron}, frameindices::Vector{<:Int}; atol = 1e-8)

Topological interlocking test by Wang. Only face face contacts are modeled.
"""
function titest(assembly::Vector{<:Polyhedron}, frameindices::Vector{<:Int}; atol = 1e-8)
    # if any block in the assembly contains flat edges, flatten the block
    for (i,block) in enumerate(assembly)
        if any([isflatedge(block, edge, atol = atol) for edge in get_edges(block)])
            assembly[i] = flattenfacets(assembly[i])
        end
    end

    contacts = contactfacets(assembly, atol = atol)

    model = Model(Ipopt.Optimizer)
    # model = Model(() -> MOA.Optimizer(Ipopt.Optimizer))
    # set_silent(model)
    # set_optimizer_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
    # set_optimizer_attribute(model, MOA.SolutionLimit(), 50)

    @variable(model, t[1:dimension(assembly[1]), 1:length(assembly)])
    @variable(model, omega[1:dimension(assembly[1]), 1:length(assembly)])

    coms = [center_of_mass(get_verts(block)) for block in assembly]
    normals = [[outward_normal(block, facet) for facet in get_facets(block)] for block in assembly]

    # constraints by paper
    for key in keys(contacts)
        i = key[1]
        j = key[2]

        if i in frameindices && j in frameindices continue end

        for k in 1:size(contacts[key])[1]
            for l in 1:size(contacts[key])[2]
                if length(contacts[key][k,l]) == 0 continue end
                
                if affinedim(contacts[key][k,l]) == 2
                    n = normals[i][k]

                    for c in contacts[key][k,l]
                        ri = c - coms[i]
                        vi = t[:,i] + cross(omega[:,i], ri)

                        rj = c - coms[j]
                        vj = t[:,j] + cross(omega[:,j], rj)

                        @constraint(model, dot(vj-vi, n) >= 0)
                    end

                elseif affinedim(contacts[key][k,l]) == 1


                else continue
                end
            end
        end
    end

    # frame constraint
    for i in frameindices
        @constraint(model, t[:,i] == 0)
        @constraint(model, omega[:,i] == 0)
    end

    @objective(model, Max, sum(t) + sum(omega))

    optimize!(model)

    # @info "t: $(display(value.(t)))"
    # @info "omega: $(display(value.(omega)))"

    if all([norm(vcat(value.(t[:,i]), value.(omega[:,i]))) < atol for i in 1:length(assembly)])
        @info "The assembly is topologically interlocking as there are no infinitesimal motions of blocks with norm < $(atol). (Maximum norm of infinitesimal motion: $(max([norm(vcat(value.(t[:,i]), value.(omega[:,i]))) for i in 1:length(assembly)]...)))"
    else
        i = findfirst([norm(vcat(value.(t[:,i]), value.(omega[:,i]))) >= atol for i in 1:length(assembly)])
        @warn "The test is inconclusive. The block $(i) has an infinitesimal motion with norm $(norm(vcat(value.(t[:,i]), value.(omega[:,i]))))"
    end

    return model
end


# function testmodel()
#     model = Model(HiGHS.Optimizer)
#     @variable(model, x[1:3] >= 0)
#     @objective(model, Max, x[3])
#     @constraint(model, dot(x, [1,0,1]) <= 1)
#     optimize!(model)
# end