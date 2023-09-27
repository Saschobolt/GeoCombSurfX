function readPolyhedronFromSTL(file::String,prints::Bool)
  input=open(file)
  temp_facets=[]
  facets=[]
  edges=[]
  verts=[]
  line=readline(input)
  while !eof(input)
    if length(split(line,"normal"))>1
      line=readline(input)
      line=readline(input)
      facet=[]
      for i in 1:3
        temp=split(line,"x ")[2]
        temp=split(temp," ")
        temp=map(g->parse(Float64,temp[g]),[1,2,3])
        push!(facet,temp)
        line=readline(input)
      end
      push!(temp_facets,facet)
      if prints 
        if (length(temp_facets) % 1000)==0 
          println(length(temp_facets))
        end
      end
    end
    line=readline(input)
  end
  push!(verts,temp_facets[1][1])
  for i in 1:length(temp_facets)
    push!(facets,[1,1,1])
  end
  for i in 1:length(facets)
    if prints 
      if (i % 1000)==0 
        println(i)
      end
    end
    for j in 1:3
      v=temp_facets[i][j]
      newVertex=true
      ind=1
      len=length(verts)
      while newVertex && ind<=len
        if norm(verts[ind]-v)<=10.0^(-7)
          newVertex=false
        end
        ind=ind+1
      end
      if newVertex
        push!(verts,v)
        facets[i][j]=ind
      else 
        facets[i][j]=ind-1
      end
    end
  end
  close(input)
  verts=map(v->Vector{Float64}([v[1],v[2],v[3]]),verts)
  facets=map(f->Vector{Int64}([f[1],f[2],f[3]]),facets)
  return PolyhedronByVerticesInFacets(verts,facets)
end

function readAssemblyFromSTL(file::String)
  p=readPolyhedronFromSTL(file)
  connect=connectedComponents(p)
  connect=map(pp->reArrangePolyhedron(pp),connect)
  println("need to define a frame")
  return connect
end


function readPolyhedronFromTextfile(file::String;prints=false)
  verts=[]
  facets=[]
  edges=[]
  input=open(file)
  line=readline(input)
  while !eof(input) 
    if length(split(line,"{"))>1 
      if length(split(line,"T"))>1 
        ##facets
        temp=split(line,"{")[2]
        temp=split(temp,"}")[1]
        temp=split(temp,";")
        push!(facets,map(g->parse(Int64,g)+1,temp))
        if prints 
          if (length(facets) %1000 )==0 
            println(length(facets))
          end
        end
      else
        ##vertices
        temp=split(line,"{")[2]
        temp=split(temp,"}")[1]
        temp=split(temp,", ")
        push!(verts,map(g->parse(Float64,g),temp))
        if prints 
          if (length(verts) %1000 )==0 
            println(length(verts))
          end
        end
      end
    end 
    line=readline(input)
  end
#  return verts,facets
  if edges==[] 
    return PolyhedronByVerticesInFacets(verts,facets)
  else
    return Polyhedron(verts,edges,facets)
  end
end


function readAssemblyFromTextfile(file::String)
  p=readPolyhedronFromTextFile(file)
  connect=connectedComponents(p)
 #connect=map(pp->reArrangePolyhedron(pp),connect)
  println("need to define a frame")
  return connect
end


#########################
###########################
function reArrangePolyhedron(p::Polyhedron)
  edges=copy(p.edges)
  facets=copy(p.facets)
  temp=vcat(p.facets...)
  temp=sort(temp)  
  verts=(p.verts)[temp]
  m=maximum(temp)

  t=Vector{Any}(undef,m)
  for i in 1:length(temp) 
    t[temp[i]]=i
  end
  ##change facets  
  for i in 1:length(facets) 
    facets[i]=map(i->t[i],facets[i])
  end

  ##change edges  
  for i in 1:length(edges) 
    edges[i]=map(i->t[i],edges[i])
  end
  return Polyhedron(verts,edges,facets)
end
