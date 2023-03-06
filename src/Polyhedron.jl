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

#########################################################################
################################## combinatorics
##########################################################################

function NumberOfVertices(poly::Polyhedron)
    return length(poly.verts)
end

function NumberOfEdges(poly::Polyhedron)
    return length(poly.edges)
end

function NumberOfFaces(poly::Polyhedron)
    return length(poly.facets)
end

function EulerCharacteristic(poly::Polyhedron)
    return NumberOfVertices(poly)-NumberOfEdges(poly)+NumberOfFaces(poly)
end


#### EdgesOfFaces

function EdgesOfFaces(poly::Polyhedron)
    res=[]
    return map(i->EdgesOfFace(poly.facets[i]),1:NumberOfFaces(poly))
end

function EdgesOfFace(poly::Polyhedron,face::Vector{Int})
    res=[]
    edges=poly.edges
    for edge in edges 
	p1=help_Position(face,edge[1])
	p2=help_Position(face,edge[2])
	if p1!= false && p2 != false 
	    if abs(p1-p2)==1 || abs(p1-p2)==length(face)-1
		push!(res,edge)
	    end
	end
    end
    return res
end

function EdgesOfFace(poly::Polyhedron,face::Int)
	return EdgesOfFace(poly,poly.facets[face])
end


####FacesOfEdges
function FacesOfEdges(poly::Polyhedron) 
    return map(i->FacesOfEdge(poly,i),1:NumberOfEdges(poly))
end

function FacesOfEdge(poly::Polyhedron,edge::Vector{Int})
    res=[]
    faces=poly.facets
    for face in faces 
        if edge[1] in face && edge[2] in face
	    p1=help_Position(face,edge[1])
	    p2=help_Position(face,edge[2])
	    if p1!= false && p2 != false 
	        if abs(p1-p2)==1 || abs(p1-p2)==length(face)-1
		    push!(res,face)
	        end
	    end
	end
    end
    return res
end

function FacesOfEdge(poly::Polyhedron, edge::Int)
	return FacesOfEdge(poly,poly.edges[edge])
end

#### FacesOfVertices
function FacesOfVertices(poly::Polyhedron )
    res=[]
    for v in 1:NumberOfVertices(poly)
	push!(res,FacesOfVertex(poly,v))
    end
    return res
end

function FacesOfVertex(poly::Polyhedron,v::Int)
    res=[]
    faces=poly.facets
    for f in faces 
        if v in f 
            push!(res,f)
        end
    end
    return res
end

####EdgesOfVertices
function EdgesOfVertex(poly::Polyhedron, v::Int)
    res=[]
    for edge in poly.edges 
	if v in edge
	    push!(res,edge)
	end
    end
    return res
end

function EdgesOfVertices(poly::Polyhedron)
    return map(i->EdgesOfVertex(poly,i),1:NumberOfVertices(poly))
end

####FaceDegreesOfVertices
function FaceDegreesOfVertices(poly::Polyhedron)
    fov=FacesOfVertices(poly)
    return map(i->length(i),fov)
end 

function FaceDegreeOfVertex(poly::Polyhedron,v::Int)
    return length(FacesOfVertex(poly,v))
end

##### Boundaery Edges
function  IsBoundaryEdge(poly::Polyhedron,edge::Vector{Int})
    foe=FacesOfEdge(poly,edge)
    return length(foe)==1 
end

function  IsBoundaryEdge(poly::Polyhedron,edge::Int)
    foe=FacesOfEdge(poly,edge)
    return length(foe)==1 
end
function BoundaryEdges(poly::Polyhedron)
    res=[]
    edges=poly.edges
    for edge in edges 
	if IsBoundaryEdge(poly,edge) 
	    push!(res, edge)
	end
    end
    return res
end

####InnerEdges

function  IsInnerEdge(poly::Polyhedron,edge::Int)
    foe=FacesOfEdge(poly,edge)
    return length(foe)==2 
end

function  IsInnerEdge(poly::Polyhedron,edge::Vector{Int})
    foe=FacesOfEdge(poly,edge)
    return length(foe)==2 
end

function InnerEdges(poly::Polyhedron)
    res=[]
    edges=poly.edges
    for edge in edges 
	if IsInnerEdge(poly,edge) 
	    push!(res, edge)
	end
    end
    return res
end


#### RamifiedEdges
function IsRamifiedEdge(poly::Polyhedron,edge::Vector{Int})
    foe=FacesOfEdge(poly,edge)
    return length(foe)>=3 
end
function IsRamifiedEdge(poly::Polyhedron,edge::Int)
    foe=FacesOfEdge(poly,edge)
    return length(foe)>=3 
end
function RamifiedEdges(poly::Polyhedron)
    res=[]
    edges=poly.edges
    for edge in edges
	if IsRamifiedEdge(poly,edge) 
	    push!(res, edge)
	end
    end
    return res
end

###IsClosedSurface
function IsClosedSurface(poly::Polyhedron)
    return poly.edges==InnerEdges(poly)
end


#### IsCOnnected ( TODO needs to be tested)
function IsConnected(poly::Polyhedron)
    faces=poly.facets
    visitedFaces=[faces[1]]
    visitedEdges=EdgesOfFace(poly,faces[1])
    faces=filter(f->f!=faces[1],faces)
    tempFaces=0
    while tempFaces != visitedFaces 
	tempF=visitedFaces
	tempE=visitedEdges
	for edge in tempE
	    foe=FacesOfEdge(poly,edge)
	    for f in foe 
		if !(f in visitedFaces)
		    push!(tempF,f)
		    for edge2 in EdgesOfFace(poly,f) 
			if !(edge2 in visitedEdges) && !([edge2[2],edge2[1]] in visitedEdges) 
			    push!(tempE,edge2)
			end
		    end
		end
	    end
	end
	tempFaces=visitedFaces
	visitedFaces=tempF
	visitedEdges=tempE
    end
    return NumberOfFaces(poly)==length(visitedFaces)
end






#### BoundaryVertex
function IsBoundaryVertex(poly::Polyhedron,v::Int)
    eov=EdgesOfVertex(poly,v)
    boundary=[]
    for edge in eov
        if IsBoundaryEdge(poly,edge) 
	    push!(boundary,edge)
	end
    end
    return boundary!=[]
end

function BoundaryVertices(poly::Polyhedron)
    filter(v->IsBoundaryVertex(poly,v),1:NumberOfVertices(poly))
end

####InnerVertex
function IsInnerVertex(poly::Polyhedron,v::Int)
    eov=EdgesOfVertex(poly,v)
    boundary=[]
    for edge in eov
        if IsBoundaryEdge(poly,edge) 
	    push!(boundary,edge)
	end
    end
    return boundary==[]
end

function InnerVertices(poly::Polyhedron)
    filter(v->IsInnerVertex(poly,v),1:NumberOfVertices(poly))
end


#TODO
#function IsRamifiedVertex(poly::Polyhedron,v::Int)
#end

#function RamifiedVertices()
#end
#TODO
#function OrientatePolyhedron(poly::Polyhedron)
#end


#TODO calling this function results in killing my terminal 
function DisjointUnion(poly::Polyhedron,poly2::Polyhedron)
    num=NumberOfVertices(poly)
    vertices=poly.verts

    for v in poly2.verts 
	push!(vertices,v)
    end

    edges=poly.edges
    for edge in poly2.edges
	push!(edges,[num+edge[1],num+edge[2]])
    end

    faces=poly.facets
    for face in poly2.facets
	push!(faces,[num+face[1],num+face[2],num+face[3]])
    end
    return (vertices,edges,faces)
end
#helperfunctions

function help_Position(l::Vector{Vector{Int}},x::Vector{Int})
    x2=[]
    for i in 1:length(x)
        push!(x2,x[length(x)-i+1])
    end
    for i in 1:length(l)
	if l[i]==x || l[i]==x2
	    return i
	end
    end
    return false
end

function help_Position(l::Vector{Vector{Int}},x::Int)
    for i in 1:length(l)
	if l[i]==x 
	    return i
	end
    end
    return false
end

function help_Position(l::Vector{Int},x::Int)
    for i in 1:length(l)
	if l[i]==x 
	    return i
	end
    end
    return false
end

function help_OrientateFace(poly::Polyhedron,face::Vector{Int},edge::Vector{Int})
    orientFace=[]
    p1=help_Position(face,edge[1])
    p2=help_Position(face,edge[2])
    print(p1," ",p2,"\n")
    if p1-p2==-1 || p1-p2==length(face)-1
        for i in 1:length(face)
	    push!(orientFace,face[length(face)-i+1])
	end 
    else
        orientFace=face
    end

    return orientFace

end

