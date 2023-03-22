#TODO change all "faces" into "facets"
function ComputeIntersections(poly1::Polyhedron,poly2::Polyhedron)
    res=[]
    for face1 in poly1.facets
        for face2 in poly1.facets
	    intersect:=[]
	    coordinates1=map(i->poly1.verts[i],face1)
	    coordinates2=map(i->poly2.verts[i],face2)
	    v=coordinates1[face1[1]]
	    v1=coordinates1[face1[2]]-coordinates1[face1[1]]
	    v2=coordinates1[face1[3]]-coordinates1[face1[1]]
	    v3=help_CrossProduct(v1,v2)
	    w=coordinates2[face2[1]]
	    w1=coordinates2[face2[2]]-coordinates2[face2[1]]
	    w2=coordinates2[face2[3]]-coordinates2[face2[1]]
	    mat=(Matrix([v1[1],v2[1],v3[1];v1[2],v2[2],v3[2];v1[3],v2[3],v3[3]]))^(-1)
	    sol=mat*(w-v)
	    if sol[3]^2<=0.0001
		temp=map(i->inpolygon3d(coordinates1[i], coordinates2,1e-8)
		temp2=map(i->inpolygon3d(coordinates2[i], coordinates1,1e-8) 		
		## face1 contained in face2 
		if temp==[1,1,1]
		    push!(res,coordinates1)
		## face 2 contained in face1 
		elseif temp2==[1,1,1]
		    push!(res,coordinates2) 
		else
			for e1 in EdgesOfFace(poly1)
			    for e2 in EdgesOfFace(poly2)
				v=coordinate[e]
			    end
			end 
		end
	    end
	end
    end 
    return res
end



