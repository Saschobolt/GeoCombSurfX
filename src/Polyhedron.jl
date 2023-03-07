using Plots
using PolygonOps


struct Polyhedron
    verts::Vector{Vector{Float32}} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{Vector{Int}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{Vector{Int}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
end


# plot Polyhedron
function plotPolyhedron!(p::Plots.Plot, poly::Polyhedron; fill_color = :red, kwargs...)
    for facet in poly.facets
        xvals = append!(map(ind -> verts[ind][1], facet), verts[facet[1]][1])
        yvals = append!(map(ind -> verts[ind][2], facet), verts[facet[1]][2])
        zvals = append!(map(ind -> verts[ind][3], facet), verts[facet[1]][3])

        plot!(p, xvals, yvals, zvals, kwargs...)
    end
    return p
end


# triangulate surface of Polyhedron
"""

"""
function triangulatePolyhedron(poly::Polyhedron)::Polyhedron
    newVerts = poly.verts
    newEdges = poly.edges
    newFacets = []


    for facet in poly.facets
        subfacet = facet

        while length(subfacet) > 3
            if inpolygon((poly.verts[subfacet[length(subfacet)]] + poly.verts[subfacet[2]]) / 2, append!(map(i -> poly.verts[i], subfacet), [poly.verts[subfacet[1]]])) == 1
                append!(newEdges, [(subfacet[length(subfacet)], subfacet[2])])
                append!(newFacets, [[subfacet[length(subfacet)], subfacet[1], subfacet[2]]])
                subfacet = subfacet[2:length(subfacet)]
                continue
            elseif inpolygon((poly.verts[subfacet[length(subfacet) - 1]] + poly.verts[subfacet[1]]) / 2, append!(map(i -> poly.verts[i], subfacet), [poly.verts[subfacet[1]]])) == 1
                append!(newEdges, [(subfacet[length(subfacet) - 1], subfacet[1])])
                append!(newFacets, [[subfacet[length(subfacet) - 1], subfacet[length(subfacet)], subfacet[1]]])
                subfacet = subfacet[1:length(subfacet) - 1]
                continue
            else
                for ind in 2:length(subfacet) - 1
                    if inpolygon((poly.verts[subfacet[ind - 1]] + poly.verts[subfacet[ind + 1]]) / 2, append!(map(i -> poly.verts[i], subfacet), [poly.verts[subfacet[1]]])) == 1
                        append!(newEdges, [(subfacet[ind - 1], subfacet[ind + 1])])
                        append!(newFacets, [[subfacet[ind - 1], subfacet[ind], subfacet[ind + 1]]])
                        subfacet = subfacet[1:end .!= ind]
                        break
                    end
                end
            end
        end

        append!(newFacets, [subfacet])
    end

    print(newEdges, "\n")
    print(newFacets)

    return Polyhedron(newVerts, newEdges, newFacets)
end

#########################################################################
################################## rigid motions
##########################################################################

#TODO
#function TranslatePolyhedron() 
#end

#function MirorPolyhedron()
#end

#function RotatePolyhedron()
#end


