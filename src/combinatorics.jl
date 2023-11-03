#include("Polyhedron.jl")

#########################################################################
#  TODO: Methode, welche die Facets von Polyeder so umschreibt, dass aufeinanderfolgende Vertizes in Facet durch Kante verbunden sind (siehe bspw. Cuboctahedron, Facet [1,3,5,7])
#########################################################################

#########################################################################
################################## combinatorics
##########################################################################

function NumberOfVertices(poly::AbstractPolyhedron)
    return length(poly.verts)
end

function NumberOfEdges(poly::AbstractPolyhedron)
    return length(poly.edges)
end

function NumberOfFacets(poly::AbstractPolyhedron)
    return length(poly.facets)
end

function EulerCharacteristic(poly::AbstractPolyhedron)
    return NumberOfVertices(poly)-NumberOfEdges(poly)+NumberOfFacets(poly)
end


#### EdgesOfFaces

function EdgesOfFacets(poly::AbstractPolyhedron)
    res=[]
    return map(i->EdgesOfFacet(poly.facets[i]),1:NumberOfFacets(poly))
end

function EdgesOfFacet(poly::AbstractPolyhedron,facet::Vector{Int})
    res=[]
    for i in 1:length(facet)-1
	push!(res,[facet[i],facet[i+1]])
    end
    push!(res,[facet[length(facet)],facet[1]])
    return res
end

function EdgesOfFacet(poly::AbstractPolyhedron,facet::Int)
	return EdgesOfFacet(poly,poly.facets[facet])
end


####FacesOfEdges
function FacetsOfEdges(poly::AbstractPolyhedron) 
    return map(i->FacetsOfEdge(poly,i),1:NumberOfEdges(poly))
end

function FacetsOfEdge(poly::AbstractPolyhedron,edge::Vector{<:Int})
    res=[]
    facets=poly.facets
    for facet in facets 
        if edge[1] in facet && edge[2] in facet
	    p1=help_Position(facet,edge[1])
	    p2=help_Position(facet,edge[2])
	    if p1!= false && p2 != false 
	        if abs(p1-p2)==1 || abs(p1-p2)==length(facet)-1
		    push!(res,facet)
	        end
	    end
	end
    end
    return res
end

function FacetsOfEdge(poly::AbstractPolyhedron, edge::Int)
	return FacetsOfEdge(poly,poly.edges[edge])
end

#### FacesOfVertices
function FacetsOfVertices(poly::AbstractPolyhedron )
    res=[]
    for v in 1:NumberOfVertices(poly)
	push!(res,FacetsOfVertex(poly,v))
    end
    return res
end

function FacetsOfVertex(poly::AbstractPolyhedron,v::Int)
    res=[]
    facets=poly.facets
    for f in facets 
        if v in f 
            push!(res,f)
        end
    end
    return res
end

####EdgesOfVertices
function EdgesOfVertex(poly::AbstractPolyhedron, v::Int)
    res=[]
    for edge in poly.edges 
	if v in edge
	    push!(res,edge)
	end
    end
    return res
end

function EdgesOfVertices(poly::AbstractPolyhedron)
    return map(i->EdgesOfVertex(poly,i),1:NumberOfVertices(poly))
end


####FacetDegreesOfVertices
function FacetDegreesOfVertices(poly::AbstractPolyhedron)
    fov=FacetsOfVertices(poly)
    return map(i->length(i),fov)
end 

function FacetDegreeOfVertex(poly::AbstractPolyhedron,v::Int)
    return length(FacetsOfVertex(poly,v))
end

##### Boundaery Edges
function  IsBoundaryEdge(poly::AbstractPolyhedron,edge::Vector{<:Int})
    foe=FacetsOfEdge(poly,edge)
    return length(foe)==1 
end

function  IsBoundaryEdge(poly::AbstractPolyhedron,edge::Int)
    foe=FacetsOfEdge(poly,edge)
    return length(foe)==1 
end
function BoundaryEdges(poly::AbstractPolyhedron)
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

function  IsInnerEdge(poly::AbstractPolyhedron,edge::Int)
    foe=FacetsOfEdge(poly,edge)
    return length(foe)==2 
end

function  IsInnerEdge(poly::AbstractPolyhedron,edge::Vector{<:Int})
    foe=FacetsOfEdge(poly,edge)
    return length(foe)==2 
end

function InnerEdges(poly::AbstractPolyhedron)
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
function IsRamifiedEdge(poly::AbstractPolyhedron,edge::Vector{<:Int})
    foe=FacetsOfEdge(poly,edge)
    return length(foe)>=3 
end
function IsRamifiedEdge(poly::AbstractPolyhedron,edge::Int)
    foe=FacetsOfEdge(poly,edge)
    return length(foe)>=3 
end
function RamifiedEdges(poly::AbstractPolyhedron)
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
function IsClosedSurface(poly::AbstractPolyhedron)
    return poly.edges==InnerEdges(poly)
end



#### IsCOnnected 
function helpunion(l::Vector)
  res=[]
  for g in l 
   union!(res,g)
  end
  return res
end

function isconnected(poly::AbstractPolyhedron)
  facets=get_facets(poly)
  vertices=helpunion(facets)
  visitedV=[vertices[1]]
  tempV=0
  while tempV != visitedV
    tempV=deepcopy(visitedV)
    temp=Set(deepcopy(visitedV))
    visitedV=filter(f->intersect!(Set(deepcopy(f)),temp) != Set([]) ,facets)
    visitedV=helpunion(visitedV) 
  end
  return length(vertices)==length(visitedV)
end

function connectedComponents(poly::AbstractPolyhedron)
  resV=[]
  resCoor=[]
  resF=[]
  tres=[]
  verts=get_verts(poly)
  facets=get_facets(poly)
  vertices=helpunion(facets)
  remVertices=deepcopy(vertices)
  while length(remVertices) != 0
    visitedV=[remVertices[1]]
    vF=0
    tempV=0
    while tempV != visitedV
      tempV=deepcopy(visitedV)
      temp=Set(deepcopy(visitedV))
      vF=filter(f->intersect!(Set(deepcopy(f)),temp) != Set([]) ,facets)
      visitedV=helpunion(vF)
    end
    vV=map(i->verts[i],visitedV)
    push!(resV,visitedV)
    push!(tres,Vector{Vector{Float64}}(undef,maximum(visitedV)))
    push!(resCoor,vV)
    push!(resF,vF)
    filter!(v->!(v in visitedV), remVertices)
  end
  for i in 1:length(resV)
    for j in resV[i] 
      tres[i][j]=verts[j]
    end
  end
  return map(i->PolyhedronByVerticesInFacets(tres[i],resF[i]),1:length(resV))
end


#### IsCOnnected 




#### BoundaryVertex
function IsBoundaryVertex(poly::AbstractPolyhedron,v::Int)
    eov=EdgesOfVertex(poly,v)
    boundary=[]
    for edge in eov
        if IsBoundaryEdge(poly,edge) 
	    push!(boundary,edge)
	end
    end
    return boundary!=[]
end

function BoundaryVertices(poly::AbstractPolyhedron)
    filter(v->IsBoundaryVertex(poly,v),1:NumberOfVertices(poly))
end

####InnerVertex
function IsInnerVertex(poly::AbstractPolyhedron,v::Int)
    eov=EdgesOfVertex(poly,v)
    boundary=[]
    for edge in eov
        if IsBoundaryEdge(poly,edge) 
	    push!(boundary,edge)
	end
    end
    return boundary==[]
end

function InnerVertices(poly::AbstractPolyhedron)
    filter(v->IsInnerVertex(poly,v),1:NumberOfVertices(poly))
end

## Orientation
function OrientatePolyhedron(poly::AbstractPolyhedron)
  facets=get_facets(poly)
  visitedFacets=[facets[1]]
  range=collect(1:length(facets))
  res=map(i->[1,1,1],range)
  res[1]=facets[1]
  visitedEdges=EdgesOfFacet(poly,facets[1])
  tempFacets=0
  while NumberOfFacets(poly) != length(visitedFacets) 
    tempF=visitedFacets
    tempE=visitedEdges
    for edge in tempE
      foe=FacetsOfEdge(poly,edge)
      for f in foe 
        if !(f in visitedFacets)
          push!(tempF, f)
          if edge in EdgesOfFacet(poly,f)
	    ff=help_OrientateFacet(poly,f,edge)
          else
            ff=f
          end
          temp=filter(i->Set(facets[i])==Set(ff),range)
          res[temp[1]]=ff
          for edge2 in EdgesOfFacet(poly,ff) 			
	    if !(edge2 in visitedEdges) && !([edge2[2],edge2[1]] in visitedEdges) 
	      push!(tempE,edge2)
	    end
          end
        end
      end
    end
    visitedFacets=tempF
    visitedEdges=tempE
  end
  return Polyhedron(poly.verts,poly.edges,res)
end

####NeighbourFaceByEdge

function NeighbourFacetByEdge(poly::AbstractPolyhedron,facet::Vector{<:Int},edge::Vector{<:Int})
   res=[]
   for f in poly.facets
        if edge[1] in f && edge[2] in f 
	    p1=help_Position(facet,edge[1])
	    p2=help_Position(facet,edge[2])
	    if (abs(p1-p2)==1 || abs(p1-p2)==length(facet)-1) && f !=facet
		push!(res,f)
	    end
	end 
   end
   return res
end

function NeighbourFacetByEdge(poly::AbstractPolyhedron,facet::Int,edge::Int)
	return NeighbourFacetByEdge(poly,poly.facets[facet],poly.edges[edge])
end

####Disjoint Union
function DisjointUnion(poly::AbstractPolyhedron,poly2::AbstractPolyhedron)
    num=NumberOfVertices(poly)
    p=CopyPolyhedron(poly)

    for v in poly2.verts
	push!(p.verts,v)
    end
    for edge in poly2.edges
	push!(p.edges,[num+edge[1],num+edge[2]])
    end
    for facet in poly2.facets
	push!(p.facets,[num+facet[1],num+facet[2],num+facet[3]])
    end
    return p
end

## face defining vertex cycles
function FacetDefiningVertexCycles(poly::AbstractPolyhedron)
  res=[]
  for facet in poly.facets 
    facetCycle=[facet[1]]
    len=1
    while len != length(facet)
      v=facetCycle[len]
      eov=EdgesOfVertex(poly,v)
      temp=filter(v2->[v,v2] in eov || [v2,v] in eov, setdiff(facet,faceCycle))
      push!(facetCycle,temp[1])
      len=len+1
    end    
    push!(res,facetCycle)
  end
  return res
end

################helperfunctions

function help_Position(l::Vector{Vector{<:Int}},x::Vector{<:Int})
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

function help_Position(l::Vector{<:Int},x::Int)
    for i in 1:length(l)
	if l[i]==x 
	    return i
	end
    end
    return false
end

function help_OrientateFacet(poly::AbstractPolyhedron,facet::Vector{<:Int},edge::Vector{<:Int})
    orientFacet=[]
    p1=help_Position(facet,edge[1])
    p2=help_Position(facet,edge[2])
    if p1-p2==-1 || p1-p2==length(facet)-1
        for i in 1:length(facet)
	    push!(orientFacet,facet[length(facet)-i+1])
	end 
    else
        orientFacet=facet
    end

    return Vector{Int}(orientFacet)

end

function CopyPolyhedron(poly::AbstractPolyhedron)
    vertices=map(i->copy(poly.verts[i]),1:NumberOfVertices(poly))
    edges=map(i->copy(poly.edges[i]),1:NumberOfEdges(poly))
    facets=map(i->copy(poly.facets[i]),1:NumberOfFacets(poly))
    return Polyhedron(vertices,edges,facets)

end
