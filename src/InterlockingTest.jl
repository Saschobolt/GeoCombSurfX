#include("Polyhedron.jl")
#include("combinatorics.jl")
#include("contactfacets.jl")

function orientatedNormals(poly::Polyhedron)
  p=OrientatePolyhedron(poly)
  facets=get_facets(p)
  vertices=get_verts(p)
  temp=map(f->[vertices[f[2]]-vertices[f[1]],
               vertices[f[3]]-vertices[f[1]]],facets)
  normals=map(t->crossProduct(t[1],t[2]),temp)
  n=(1/norm(normals[1]))*normals[1]
  w=1/length(facets[1])*sum(map(i->vertices[i],facets[1]))
  tpoints=[]
  for i in 2:length(facets)
    f=facets[i]
    v=vertices[f[2]]
    vec2=vertices[f[1]]-v
    vec3=vertices[f[3]]-v
    mat=Matrix([vec2 vec3 (-1)*n])
    if det(mat)^2>=1e-8 
      mat=inv(mat)
      sol=mat*(w-v)
      res=w+sol[3]*n
      temp=deepcopy(f)
      push!(temp,f[1])
      if sol[3]>=1e-8 && inpolygon3d(res, map(i->vertices[i],temp)) !=0
        push!(tpoints,res)
      end
    end
  end ##
  points=[]
  for v in tpoints 
    if filter(w-> norm(w-v)<=1e-8,points) ==[] 
      push!(points,v)
    end
  end
  if length(points) % 2 != 0
    normals=map(v->(-1)*v,normals)
  end
  return normals
end


function InterlockingTestByWang2(poly::Polyhedron,frame::Vector{Int})
  # construct lp and add the variables
  assembly=connectedComponents(poly)
  indices=filter(i->!(i in frame),1:length(assembly))
  model=Model(HiGHS.Optimizer)
  @objective(model,Max,0)
  @variable(model,omega[i=indices, j=1:3])
  @variable(model,t[i=indices, j=1:3])
  println(model)
  ##edge-edge contact 
  edge_normals=[]
  for ii in 1:length(indices)-1
    for jj in ii+1:length(indices)
      i=indices[ii]
      j=indices[jj]
      p1=assembly[i]
      p2=assembly[j]
      facets1=get_facets(p1)
      facets2=get_facets(p2)
      verts1=get_verts(p1)
      verts2=get_verts(p2)
      vertices1=map(l->p1.verts[l],helpunion(facets1))
      vertices2=map(l->p2.verts[l],helpunion(facets2))
      edges1=get_edges(p1)
      edges2=get_edges(p2)
      center1=1/length(vertices)*sum(vertices1)
      center2=1/length(vertices2)*sum(vertices2)
      for e1 in p1.edges
        for e2 in p2.edges
          v=verts1[e1[1]]
          vec=verts1[e1[2]]-v
          w=verts2[e2[1]]
          vec2=verts2[e2[2]]-w
          n=crossProduct(vec,vec2)
          if norm(n)>=1e-8
            mat=transpose(Matrix([vec vec2 (-1)*n]))
            if det(mat)^2>=1e-8
              mat=inv(mat)
              sol=mat*(w-v)
              if sol[1]>=-1e-8 && sol[1]<=1.0+ 1e-8 && 
                 sol[2]>=-1e-8 && sol[2]<=1.0+ 1e-8 && 
                 sol[3]<=1e-8 && sol[3]>=-1e-8 
                c=v+sol[1]*vec   ## notation from paper 
                rc_1=c-center1
                rc_2=c-center2
                vc_1=[t[i,1]+omega[i,2]*rc_1[3]-omega[i,3]*rc_1[2],
                      t[i,2]+omega[i,3]*rc_1[1]-omega[i,1]*rc_1[3],
                      t[i,3]+omega[i,1]*rc_1[2]-omega[i,2]*rc_1[1]]
                vc_2=[t[j,1]+omega[j,2]*rc_2[3]-omega[j,3]*rc_2[2],
                      t[j,2]+omega[j,3]*rc_2[1]-omega[j,1]*rc_2[3],
                      t[j,3]+omega[j,1]*rc_2[2]-omega[j,2]*rc_2[1]]
                if transpose(((c-center1)-(c-center2)))*n>0.
                 #add equation to the lp
                  @constraint(model,(vc_1[1]-vc_2[1])*n[1]+
                                   (vc_1[2]-vc_2[2])*n[2]+
                                   (vc_1[3]-vc_2[3])*n[3]>=0)
                else
                  # add equation to the lp
                  @constraint(model,(vc_1[1]-vc_2[1])*n[1]+
                                   (vc_1[2]-vc_2[2])*n[2]+
                                   (vc_1[3]-vc_2[3])*n[3]<=0)
                end
              end
            end
          end
        end
      end
    end
  end
  ##face-face contact 
  for i in 1:length(assembly)-1
    for j in i+1:length(assembly)
      contactfacets=contactFacets2(assembly[i],assembly[j])
      if contactfacets != []
        p1=OrientatePolyhedron(assembly[i])
        p2=OrientatePolyhedron(assembly[j])
        verts1=get_verts(p1)
        verts2=get_verts(p2)
        vertices1=map(l->p1.verts[l],helpunion(facets1))
        vertices2=map(l->p2.verts[l],helpunion(facets2))
        center1=1/length(vertices1)*sum(vertices1)
        center2=1/length(vertices2)*sum(vertices2)
        temp=map(f->[vertices[f[2]]-vertices[f[1]],
                     vertices[f[3]]-vertices[f[1]]],block1.facets)
        normals1=map(t->crossProduct(t[1],t[2]),temp)
        for f in contactfacets 
          for c in f  
            rc_1=c-center1
            rc_2=c-center2
            vc_1=[t[i,1]+omega[i,2]*rc_1[3]-omega[i,3]*rc_1[2],
                  t[i,2]+omega[i,3]*rc_1[1]-omega[i,1]*rc_1[3],
                  t[i,3]+omega[i,1]*rc_1[2]-omega[i,2]*rc_1[1]]
            vc_2=[t[j,1]+omega[j,2]*rc_2[3]-omega[j,3]*rc_2[2],
                  t[j,2]+omega[j,3]*rc_2[1]-omega[j,1]*rc_2[3],
                  t[j,3]+omega[j,1]*rc_2[2]-omega[j,2]*rc_2[1]]
            if transpose(((point-center1)-(point-center2)))*n>0.
             #add equation to the lp
              @constraint(model,(vc_1[1]-vc_2[1])*n[1]+
                               (vc_1[2]-vc_2[2])*n[2]+
                               (vc_1[3]-vc_2[3])*n[3]>=0)
            else
              # add equation to the lp
              @constraint(model,(vc_1[1]-vc_2[1])*n[1]+
                               (vc_1[2]-vc_2[2])*n[2]+
                               (vc_1[3]-vc_2[3])*n[3]<=0)
            end
          end
        end
      end
    end
  end

  optimize!(model)
end



#### assemblies of candyblocks

function assembly(n1::Integer, n2::Integer, nmerges::Integer,assemblyType::Integer,k::Integer) 
  if n2 !=3 && assemblyType > 2 
    return "assembly not possible"
  end
  alpha=2*3.14159/n1
  mat=Matrix([[cos(alpha),sin(alpha),0] [-sin(alpha),cos(alpha),0] [0,0,1]])
  mat=transpose(mat)
  mat2=mat^(n1-1)
  p=tiblock(n1,n2,nmerges)
  vof=get_facets(p)
  max=length(helpunion(p.facets))
  coor=get_verts(p)
  verticesOfFacets=Vector{Vector{Int}}(undef,0)
  coordinates=Vector{Vector{Float64}}(undef,0)
  if assemblyType==1
    for i in 1:k 
      append!(verticesOfFacets,map(g->max*(i-1)*map(h->1,1:length(g))+g,vof))
      if i % 2 == 1  
        append!(coordinates,map(g->g+[0.,0.,(i-1.)],coor))
      else 
        append!(coordinates,map(g->mat*g+[0.,0.,(i-1.)],coor))
      end
    end
  elseif assemblyType==2  
    for i in 1:k 
  append!(verticesOfFacets,map(g->max*(i-1)*map(h->1,1:length(g))+g,vof))
      if i % 3==1 
        append!(coordinates,map(g->g+(i-1)*[0,0,1],coor))
      elseif i % 3==2 
        append!(coordinates,map(g->mat*g+(i-1)*[0,0,1],coor))
      else 
        append!(coordinates,map(g->mat2*g+(i-1)*[0,0,1],coor))
      end
    end
  elseif assemblyType==3  
    for i in 1:k 
      append!(verticesOfFacets,map(g->max*(i-1)*map(h->1,1:length(g))+g,vof))
      if i % 3== 1 
        append!(coordinates,map(g->g+(i-1)*[0,0,1],coor))
      elseif i % 3 ==2 
        append!(coordinates,map(g->mat*g+(i-1)*[0,0,1],coor))
      else 
        append!(coordinates,map(g->mat2*g+(i-1)*[0,0,1],coor))
      end
    end
    vec=coor[1]+1/2*(coor[2]-coor[1])
    norm=sqrt(vec[1]^2+vec[2]^2)
    vec=(2*norm+sqrt(3.)/2)/norm*vec
    for i in 1:(k-1)/3 
      append!(verticesOfFacets,map(g->(max*k+max*(i-1))*map(h->1,1:length(g))+g,vof))
      append!(coordinates,map(g->g+(i-1)*[0,0,3]+[0,0,1.5]+vec,coor))
    end
  elseif assemblyType==4 
    vec=coor[1]+1/2*(coor[2]-coor[1])
    norm=sqrt(vec[1]^2+vec[2]^2)
    vec=1/norm*vec
    a=sqrt(3.)+4*norm
    for l in 1:3 
      for i in 1:k 
	a2=((l-1)*k+(i-1))*max
        append!(verticesOfFacets,map(g->a2*map(h->1,1:length(g))+g,vof))
        if i % 3== 1  
          append!(coordinates,map(g->g+(i-1)*[0.,0.,1.]+(l-1)*a*vec,coor))
        elseif i % 3 ==2
          append!(coordinates,map(g->mat*g+(i-1)*[0.,0.,1.]+(l-1)*a*vec,coor))
        else 
          append!(coordinates,map(g->mat2*g+(i-1)*[0.,0.,1.]+(l-1)*a*vec,coor))
        end
      end
    end
    a2=(2*norm+sqrt(3.)/2)
    for l in 1:2 
      for i in 1:(k-1)/3 
        append!(verticesOfFacets,
          map(g->length(union!(verticesOfFacets))*map(h->1,1:length(g))+g,vof))
        append!(coordinates,
          map(g->g+(i-1)*[0.,0.,3.]+[0.,0.,1.5]+((l-1)*a+a2)*vec,coor))
      end
    end
  end
  return PolyhedronByVerticesInFacets(coordinates,verticesOfFacets)
end


function PolyhedronByVerticesInFacets(verts::Vector,facets::Vector)
  edges=Vector{Vector{Int}}(undef,0)
  for f in facets
    temp=map(i->[f[i],f[i+1]],1:length(f)-1)
    push!(temp,[f[1],f[length(f)]])
    for ee in temp
      if filter(e->Set(e)==Set(ee) ,edges)==[]
        push!(edges,ee)
      end
    end
  end
  return Polyhedron(verts,edges,facets)
end



function cubeAssembly(n::Integer)
  verts=get_verts(Cube)
  vertices=Vector{Vector{Float64}}(undef,0)
  fac=get_facets(Cube)
  facets=Vector{Vector{Int}}(undef,0)
  println("fac",fac)
  for i in 1:n
    for j in 1:n
      append!(vertices,map(g->[i-1,j-1,0]+g,verts))
    end
  end
  for i in 1:n*n
    append!(facets,map(f->8*(i-1)*[1,1,1,1]+f,fac))
  end
  return PolyhedronByVerticesInFacets(vertices,facets)
end


function contactFacets2(p1::Polyhedron,p2::Polyhedron)
  res=[]
  facets1=get_facets(p1)
  facets2=get_facets(p2)
  vertices1=get_verts(p1)
  vertices2=get_verts(p2)
  for f1 in facets1  
    for f2 in facets2 
      cf=[]
      f=[f1,f2]
      verts=[vertices1,vertices2]
      normal=0
      for i in 1:2 
        temp=deepcopy(f[3-i])
        push!(temp,f[3-i][1])
        temp=map(j->verts[3-i][j],temp)
        for j in f[i]
          v=verts[i][j]
          if inpolygon3d(v,temp) != 0
            push!(cf,1.0*v)
          end 
        end
      end
      ## vertices on edges
      tempf1=deepcopy(f1)
      push!(tempf1,f1[1])
      for i in 1:length(tempf1)-1
        v=1.0*vertices1[tempf1[i]]
        vec=1.0*vertices1[tempf1[i+1]]-v
        temp=deepcopy(f2)
        push!(temp,f2[1])
        for j in 1:length(temp)-1
          w=1.0*vertices2[temp[j]]
          vec2=1.0*vertices2[temp[j+1]]-w
          n=crossProduct(vec,vec2)
          if norm(n)>=1e-8
            normal=deepcopy(n)
            mat=Matrix([vec (-1.0)*vec2 n])
            if det(mat)^2>=1e-8
              mat=inv(mat)
              sol=mat*(w-v)
              if sol[1]>=-1e-6 && sol[1]<=1.0+1e-6 &&
                 sol[2]>=-1e-6 && sol[2]<=1.0+1e-6 &&
                 sol[3]<=1e-6 && sol[3]>=-1e-6
                 push!(cf,1.0*v+sol[1]*vec)
              end  
            end
          end
        end
      end
      temp=[]
      for v in cf 
        if filter(vv->norm(vv-v)<=1e-8,temp) ==[]
          push!(temp,v)
        end
      end
      if length(temp)>=3 
        push!(res,[temp,normal])
      end
    end
  end
  return res
end


function testComputeIntersection2(num::Int)
  if num==1 
    surf1=Polyhedron([[0, 0, 0], [1, 0, 0], [0, 1, 0]],
                     [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
    surf2=Polyhedron([[0.1, 0.1, 0.0], [0.1, 0.2, 0.0], [0.2, 0.1, 0.0]],
                     [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
  elseif num==2 
    surf1=Polyhedron([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
                     [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
    surf2=Polyhedron([[0.0, 0.5, 0.0], [1.0, 0.5, 0.0], [0.5, -1.0, 0.0]],
                     [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
  elseif num==3
    surf1=Polyhedron([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.8, 0.5, 0.0]],
                     [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
    surf2=Polyhedron([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.2, 0.5, 0.0]],
                     [[1, 2], [2, 3], [1, 3]], [[1, 2, 3]])
  elseif num==4 
    surf1=Polyhedron([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
                     [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]],
                     [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
    surf2=Polyhedron([[1, 1, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
                     [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]], 
                     [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]])
  end
  return contactFacets2(surf1,surf2)
end


