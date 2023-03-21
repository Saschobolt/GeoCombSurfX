#TODO change all "faces" into "facets"
function ComputeIntersections(poly1::Polyhedron,poly2::Polyhedron)
    res=[]
    for face1 in poly1.facets
        for face2 in poly1.facets
	    intersect:=[]
	    coordinates1=map(i->poly1.verts[i],face1)
	    coordinates2=map(i->poly2.verts[i],face2)
	    v=coordinates1[1]
	    v1=coordinates1[2]-coordinates1[1]
	    v2=coordinates1[3]-coordinates1[1]
	    v3=help_CrossProduct(v1,v2)
	    w=coordinates2[1]
	    w1=coordinates2[2]-coordinates2[1]
	    w2=coordinates2[3]-coordinates2[1]
	    mat=(Matrix([v1[1],v2[1],v3[1];v1[2],v2[2],v3[2];v1[3],v2[3],v3[3]]))^(-1)
	    sol=mat*(w-v)
	    if sol[3]^2<=0.0001
		temp=map(i->inpolygon3d(coordinates1[i], coordinates2,1e-8)
		temp2=map(i->inpolygon3d(coordinates2[i], coordinates1,1e-8) 		

		if temp==[1,1,1]
		    push!(res,coordinates1)
		elseif temp2==[1,1,1]
		    push!(res,coordinates2)
		elseif temp==[-1,-1,-1]
		    push!(res,coordinates1)
		else 
##TODO	    
		elseif length(filter(i->i=1,temp))=1
			pos=help_Position(temp,1)
			v=coordinates[pos]
			temp=filter(i->i !=pos,[1,2,3])
			sol1=help_coputeIntersectionFaceEdge([v,coordinates[temp[1]]],coordintes2)
			sol2=help_coputeIntersectionFaceEdge([v,coordinates[temp[1]]],coordinates2)

		elseif length(filter(i->i=1,temp2))=1
#TODO copy from above case

		elseif length(filter(i->i=1,temp))=2
		elseif length(filter(i->i=1,temp2))=2

# TODO copy from above case 
		else
 
		end
	    end
	end
    end 
    return res
end



