using JuMP
using HiGHS
using Ipopt
using Polyhedra

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
    contacts = contactfacets(assembly, atol = atol)

    model = Model(Ipopt.Optimizer)
    @variable(model, t[1:dimension(assembly[1]), 1:length(assembly)])
    @variable(model, omega[1:dimension(assembly[1]), 1:length(assembly)])

    coms = [center_of_mass(get_verts(block)) for block in assembly]
    normals = [[outward_normal(block, facet) for facet in get_facets(block)] for block in assembly]

    # constraints by paper
    for key in keys(contacts)
        i = key[1]
        j = key[2]

        for k in 1:size(contacts[key])[1]
            for l in 1:size(contacts[key])[2]
                test = (contacts[key][k,l])

                if length(contacts[key][k,l]) == 0 continue end
                if affinedim(contacts[key][k,l]) < 2 continue end
                @info "key: $(key)"
                @info "facet: $(get_facets(assembly[i])[k])"

                n = normals[i][k]

                for c in contacts[key][k,l]
                    ri = c - coms[i]
                    vi = t[:,i] + cross(omega[:,i], ri)

                    rj = c - coms[j]
                    vj = t[:,j] + cross(omega[:,j], rj)

                    @constraint(model, dot(vj-vi, n) >= 0)
                end
            end
        end
    end

    # frame constraint
    for i in frameindices
        @constraint(model, t[:,i] == 0)
        @constraint(model, omega[:,i] == 0)
    end

    @objective(model, Max, sum(t.^2) + sum(omega.^2))

    optimize!(model)

    return model
end


# function testmodel()
#     model = Model(HiGHS.Optimizer)
#     @variable(model, x[1:3] >= 0)
#     @objective(model, Max, x[3])
#     @constraint(model, dot(x, [1,0,1]) <= 1)
#     optimize!(model)
# end