using JuMP
import HiGHS
import Polyhedra

include("../Polyhedron.jl")
include("../affine_geometry.jl")
include("../decomposition.jl")


function contacts(poly1::AbstractPolyhedron, poly2::AbstractPolyhedron; atol = 1e-8)
    coords_facets1 = [get_verts(poly1)[:, facet] for facet in get_facets(poly1)]
    coords_facets2 = [get_verts(poly2)[:, facet] for facet in get_facets(poly2)]

    poly_facets1 = [Polyhedra.polyhedron(Polyhedra.vrep(transpose(coords))) for coords in coords_facets1]
    poly_facets2 = [Polyhedra.polyhedron(Polyhedra.vrep(transpose(coords))) for coords in coords_facets2]

    intersections = transpose.(Matrix{Matrix{<:Real}}([Polyhedra.vrep(Polyhedra.intersect(face1, face2)).V for face1 in poly_facets1, face2 in poly_facets2]))

    intersections = map(entry -> map(coor -> Base.round.(coor, digits = -Int(log10(atol))), entry), intersections)
    intersections = map(entry -> map(coor -> map(x -> x == -0.0 ? 0.0 : x, coor), entry), intersections)

    return intersections
end


function contacts(assembly::Vector{<:AbstractPolyhedron}; atol = 1e-8)
    contacts_dict = Dict()

    for i in 1:length(assembly)
        for j in i+1:length(assembly)
            intersections = contacts(assembly[i], assembly[j], atol = atol)
            if any(length.(intersections) .> 0)
                get!(contacts_dict, [i,j], intersections)
            end
        end
    end 

    return contacts_dict
end

function wangtest(assembly::Vector{<:AbstractPolyhedron}, frameindices::Vector{<:Int}; basis_translation::Matrix{<:Real} = [1 0 0; 0 1 0; 0 0 1], basis_rotation::Matrix{<:Real} = [1 0 0; 0 1 0; 0 0 1], atol = 1e-8)
    @assert rank(basis_translation) == min(size(basis_translation)...) || all([i == 0 for i in basis_translation]) "The translation basis must be of full rank or zero matrix."
    @assert rank(basis_rotation) == min(size(basis_rotation)...) || all([i == 0 for i in basis_rotation]) "The rotation basis must be of full rank or zero matrix."
    
    # if any block in the assembly contains flat edges, flatten the block
    for (i,block) in enumerate(assembly)
        if any([isflatedge(block, edge, atol = atol) for edge in get_edges(block)])
            assembly[i] = flattenfacets(block, atol = atol)
        end

        assembly[i] = orient_facets(block, atol = atol)
    end

    contacts_dict = contacts(assembly, atol = 1e-12)

    model = Model(HiGHS.Optimizer)
    # model = Model(() -> MOA.Optimizer(Ipopt.Optimizer))
    # set_silent(model)
    # set_optimizer_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
    # set_optimizer_attribute(model, MOA.SolutionLimit(), 50)

    @variable(model, t[1:size(basis_translation)[2], 1:length(assembly)])
    if all([i == 0 for i in basis_translation])
        @constraint(model, t .== 0)
    end
    @variable(model, omega[1:size(basis_rotation)[2], 1:length(assembly)])
    if all([i == 0 for i in basis_rotation])
        @constraint(model, omega .== 0)
    end

    coms = [center_of_mass(get_verts(block)) for block in assembly]
    normals = [[outward_normal(block, facet) for facet in get_facets(block)] for block in assembly]

    # constraint that the velocity of c relative to block i points in the same direction as n.
    function nonpen_constraint!(i,j,n,c)
        ri = c - coms[i]
        vi = basis_translation * t[:,i] + cross(basis_rotation * omega[:,i], ri)

        rj = c - coms[j]
        vj = basis_translation * t[:,j] + cross(basis_rotation * omega[:,j], rj)

        @constraint(model, dot(vj-vi, n) >= 0)
    end
    
    # constraint for twodimensional contacts
    function twodim_constraint!(i,j,polygon)
        index = findfirst(map(entry -> size(entry)[2] == 0 ? false : Set([entry[:,col] for col in 1:size(entry)[2]]) == Set([polygon[:,col] for col in 1:size(polygon)[2]]), contacts_dict[[i,j]]))
        k = index[1]

        n = normals[i][k]
        for c in [polygon[:, col] for col in 1:size(polygon)[2]]
            nonpen_constraint!(i,j,n,c)
        end
    end


    function extaffdim(a)
        if a == []
            return 0
        end

        return affinedim(a)
    end

    # constraints by paper
    for key in keys(contacts_dict)
        # @info "key: $(key)"
        i = key[1]
        j = key[2]

        if i in frameindices && j in frameindices continue end
        
        twodim = filter(x -> extaffdim(x) == 2, contacts_dict[[i,j]])
        for polygon in twodim
            twodim_constraint!(i,j,polygon)
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
        val, i = findmax(i -> norm(vcat(value.(t[:,i]), value.(omega[:,i]))), 1:length(assembly))
        @warn "The test is inconclusive. The block $(i) has an infinitesimal motion with norm $(val)"
    end

    return model
end

"""
    titest(assembly::Vector{<:AbstractPolyhedron}, frameindices::Vector{<:Int}; atol = 1e-8)

Topological interlocking test by Wang. Only face face contacts are modeled.
"""
function titest(assembly::Vector{<:AbstractPolyhedron}, frameindices::Vector{<:Int}; atol = 1e-8)
    # if any block in the assembly contains flat edges, flatten the block
    for (i,block) in enumerate(assembly)
        if any([isflatedge(block, edge, atol = atol) for edge in get_edges(block)])
            assembly[i] = flattenfacets(block, atol = atol)
        end

        assembly[i] = orient_facets(block)
    end

    contacts_dict = contacts(assembly, atol = 1e-9)

    model = Model(HiGHS.Optimizer)
    # model = Model(() -> MOA.Optimizer(Ipopt.Optimizer))
    # set_silent(model)
    # set_optimizer_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
    # set_optimizer_attribute(model, MOA.SolutionLimit(), 50)

    @variable(model, t[1:dimension(assembly[1]), 1:length(assembly)])
    @variable(model, omega[1:dimension(assembly[1]), 1:length(assembly)])

    coms = [center_of_mass(get_verts(block)) for block in assembly]
    normals = [[outward_normal(block, facet) for facet in get_facets(block)] for block in assembly]

    # constraint that the velocity of c relative to block i points in the same direction as n.
    function nonpen_constraint!(i,j,n,c)
        ri = c - coms[i]
        vi = t[:,i] + cross(omega[:,i], ri)

        rj = c - coms[j]
        vj = t[:,j] + cross(omega[:,j], rj)

        @constraint(model, dot(vj-vi, n) >= 0)
    end
    
    # constraint for twodimensional contacts
    function twodim_constraint!(i,j,polygon)
        index = findfirst(map(entry -> size(entry)[2] == 0 ? false : Set([entry[:,col] for col in 1:size(entry)[2]]) == Set([polygon[:,col] for col in 1:size(polygon)[2]]), contacts_dict[[i,j]]))
        k = index[1]

        n = normals[i][k]
        for c in [polygon[:, col] for col in 1:size(polygon)[2]]
            nonpen_constraint!(i,j,n,c)
        end
    end
    
    # constraint for onedimensional contacts
    function onedim_constraint!(i,j,line)
        # @info "contact line: $(line)"
        # indices of contact line segment in contacts[key]. Number of those is 2 if the contact is edge-face, 4 if it is edge-edge
        indices = findall(entry -> typeof(size(entry)) <: Tuple{Integer} ? false : Set([entry[:, col] for col in 1:size(entry)[2]]) == Set([line[:,col] for col in 1:size(line)[2]]), contacts_dict[[i,j]])
        # @info "i: $(i), j: $(j), indices: $(indices)"
        # @info "indices: $(indices)"

        if length(indices) == 3
            # it can happen that two of the 4 adjacent facets of the line segment align (as the contact between those 2 facets is 2d)
            # add the contacts to the indices in that case.
            facesi = unique([indices[1][1], indices[2][1], indices[3][1]])
            facesj = unique([indices[1][2], indices[2][2], indices[3][2]])
            indices = CartesianIndex.([(facesi[1], facesj[1]), (facesi[1], facesj[2]), (facesi[2], facesj[1]), (facesi[2], facesj[2])])
        end
        
        if length(indices) == 2
            # @info "pure edge face contact: indices $(indices)"
            # pure edge face contact
            # determine which block contains the face
            if indices[1][1] == indices[2][1]
                # block i contains the face
                n = normals[i][indices[1][1]]
                for c in [line[:, col] for col in 1:size(line)[2]]
                    nonpen_constraint!(i,j,n,c)
                end

            else 
                # block j contains the face
                n = normals[j][indices[1][2]]
                for c in [line[:, col] for col in 1:size(line)[2]]
                    nonpen_constraint!(j,i,n,c)
                end
            end

        elseif length(indices) == 4
            # pure edge edge contact
            # determine edges
            facetsi_inds = unique([indices[k][1] for k in eachindex(indices)])
            facetsi = get_facets(assembly[i])[facetsi_inds]
            facetsj_inds = unique([indices[k][2] for k in eachindex(indices)])
            facetsj = get_facets(assembly[j])[facetsj_inds]
            # @info "facetsi: $(facetsi)"

            edgei = Base.intersect(incedges(assembly[i], facetsi[1]), incedges(assembly[i], facetsi[2]))[1]
            edgej = Base.intersect(incedges(assembly[j], facetsj[1]), incedges(assembly[j], facetsj[2]))[1]

            edgetypei = edgetype(assembly[i], edgei, is_oriented = true)
            edgetypej = edgetype(assembly[j], edgej, is_oriented = true)

            if edgetypei == "convex" && edgetypej == "convex"
                # @info "edge $(edgei) (convex) of block $(i) and edge $(edgej) (convex) of block $(j)"
                # TODO: Implementieren
                return
            elseif edgetypei == "concave" && edgetypej == "convex"
                # @info "edge $(edgei) (concave) of block $(i) and edge $(edgej) (convex) of block $(j)"
                n1 = normals[i][facetsi_inds[1]]
                n2 = normals[i][facetsi_inds[2]]
                for c in [line[:, col] for col in 1:size(line)[2]]
                    nonpen_constraint!(i,j,n1,c)
                    nonpen_constraint!(i,j,n2,c)
                end
            elseif edgetypei == "convex" && edgetypej == "concave"
                # @info "edge $(edgej) (concave) of block $(j) and edge $(edgei) (convex) of block $(i)"
                n1 = normals[j][facetsj_inds[1]]
                n2 = normals[j][facetsj_inds[2]]
                for c in [line[:, col] for col in 1:size(line)[2]]
                    nonpen_constraint!(j,i,n1,c)
                    nonpen_constraint!(j,i,n2,c)
                end
            end
        end
    end


    function extaffdim(a)
        if a == []
            return 0
        end

        return affinedim(a)
    end

    # constraints by paper
    for key in keys(contacts_dict)
        # @info "key: $(key)"
        i = key[1]
        j = key[2]

        if i in frameindices && j in frameindices continue end
        
        twodim = filter(x -> extaffdim(x) == 2, contacts_dict[[i,j]])
        for polygon in twodim
            twodim_constraint!(i,j,polygon)
        end
        
        onedim = filter(x -> extaffdim(x) == 1, contacts_dict[[i,j]])
        onedim = collect.(unique(Set([[line[:,col] for col in 1:size(line)[2]] for line in onedim])))
        onedim = [hcat(line...) for line in onedim]
        # @info "number of 1d contacts: $(length(onedim))"
        # @info "onedim: $(onedim)"
        for line in onedim
            onedim_constraint!(i,j,line)
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
        val, i = findmax(i -> norm(vcat(value.(t[:,i]), value.(omega[:,i]))), 1:length(assembly))
        @warn "The test is inconclusive. The block $(i) has an infinitesimal motion with norm $(val)"
    end

    return model
end