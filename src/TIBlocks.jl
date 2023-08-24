function Candy2(n::Int)
  l=Vector{Vector{Vector{Int}}}(undef,0)
  push!(l,[[1,n+1,n+2,2],[2,1,4,5]])
  push!(l,[[1,2,2*n+2,2*n+1],[2,5,4,1]]) 
  push!(l,[[n+2,n+1,2*n+1,2*n+2],[2,5,4,1]])
  push!(l,[[n+n/2+2,n/2+2,n/2+1,n+n/2+1],[2,1,4,5]])
  push!(l,[[n/2+1,n/2+2,2*n+6+1,2*n+6+2],[2,5,4,1]])
  push!(l,[[n+n/2+1,n+n/2+2,2*n+6+1,2*n+6+2],[5,2,1,4]])
  temp=Prisma(n)
  prism3=Prisma(3)
  for g in l 
    println(g[1],typeof(g[1]))
    temp=merge(temp,prism3,[g[1]],[g[2]])
  end
  return temp
end

function CubeCandy(n::Int)
  quad=Cube
  l=[[[3,4,8,7],[3,4,8,7]]]
  push!(l,[[1,2,6,5],[1,2,6,5]])
  for g in l 
    quad=merge(quad,Cube,[g[1]],[g[2]])  
  end
  l=[[[1,n+1,n+2,2,],[4,1,2,3]]]
  push!(l,[[n+n/2+2,n/2+2,n/2+1,n+n/2+1],[4,1,2,3]])
  temp=Prism(n)
  for g in l  
    temp=merging(temp,quad,[g[1]],[g[2]])
  end
  return temp
end

#für zweier
function generalizedCandy(n::Int)
  l=[[[1,n+1,n+2,2],[2,1,4,5]]]
  temp=Prisma(n)
  for i in 1:n/2
    l=[[[2*(i-1)+1,n+2*(i-1)+1,n+2*(i-1)+2,2*(i-1)+2],[2,1,4,5]]]
    push!(l,[[2*(i-1)+1,2*(i-1)+2,2*n+2*(i-1)+2,2*n+2*(i-1)+1],[2,5,4,1]]) 
    push!(l,[[n+2*(i-1)+2,n+2*(i-1)+1,2*n+2*(i-1)+1,2*n+2*(i-1)+2],[2,5,4,1]])   
    for g in l 
      temp=helpAttach(p,pr,trianPrisma,prTrianPrisma,g,-1)

    end
  end
  return temp
end


#für zweier
function generalizedCubeCandy(n::Int)	
  quad=Cube
  l=[[[3,4,8,7],[3,4,8,7]]]
  push!(l,[[1,2,6,5],[1,2,6,5]])
  copyCube=Cube
  for g in l 
    quad=merging(copyCube,quad,[g[1]],[g[2]])  
  end
  temp=Prism(n)
  for i in 1:n/2 
    g=[[2*(i-1)+1,n+2*(i-1)+1,n+2*(i-1)+2,2*(i-1)+2],[4,1,2,3]]
    temp=merging(temp,quad,[g[1]],[g[2]])
  end
  return temp
end


