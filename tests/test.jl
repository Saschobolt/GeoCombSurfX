using Test
# using .GeoCombSurfX

include("../src/GeoCombSurfX.jl")
# using .GeoCombSurfX

function test()
    @testset  "affine_geometry.jl" begin
        # "test validity of matrix functions -> rank, colspace, ..."
        A = rand(Float64, (10,20))
        @test length(indcols_indices(A)) == rank(A)
        @test size(indcols(A)) == (10, rank(A))
        @test size(colspace(A)) == size(indcols(A))
        @test length(affinebasis_indices(A)) == rank(A) + 1
        @test size(affinebasis(A)) == (size(A)[1], rank(A) + 1)
        @test affinedim(A) == rank(A)

        #test that affinemap and rigidmap work properly  # TODO: Errors testen
        preim = affinebasis(rand(Float64, (10,20)))
        im = rand(Float64, (5, size(preim)[2]))
        @test norm(affinemap(preim, im)(preim) - im) <1e-10

        preim = rand(Float64, (3, 4))
        im = preim + ones((3,4))
        @test norm(rigidmap(preim, im)(preim) - im) < 1e-10

        # test Plane and Ray
        plane = Plane([[0,0,0], [1,0,0], [0,1,0]])
        @test normalvec(plane) == [0,0,1] || normalvec(plane) == [0,0,-1]
        ray = Ray([0,0,0], [0,0,1])
        @test intersect(ray, plane) == [0,0,0]
    end

    @testset "Framework.jl" begin
        @test_throws AssertionError Graph([1,2,3], [[1,2,3]])
        @test_throws AssertionError Graph([1,2,3], [[1,2], [1,2]])
        @test_throws AssertionError Graph([1,2,3], [[1,4]])
        @test_throws AssertionError Graph([1,2,4], [[1,2]])
        @test_throws "consist of vertices of the graph" set_edges!(Graph([1,2,3], [[1,2]]), [[1,4]])

        function graph_with_components(n)
            verts = collect(1:3*n)
            edges = vcat([[[(i-1)*3 + 1, (i-1)*3 + 2], [(i-1)*3 + 2, (i-1)*3 + 3]] for i in 1:n]...)
            return Graph(verts, edges)
        end

        for n in 2:10
            @test no_concomponents(graph_with_components(n)) == n
            @test_throws "not connected" Framework(graph_with_components(n))
        end

        @test_throws "coordinate assigned to it." Framework(ones(Float64, (3,2)), [[1,3], [1,2]])
    end

    @testset "polygonal_geometry.jl" begin
        triang = [0 1 0; 0 0 1; 0 0 0]
        polygon = [0 2 0; 0 0 2; 0 0 0]
        @test intriang3d(triang, [0,0,0]) == -1
        @test intriang3d(triang, [1/2, 0, 0]) == -1
        @test intriang3d(triang, center_of_mass(triang)) == 1
        @test intriang3d(triang, [1,1,1]) == 0
        @test intriang3d(polygon, [0.5, 0.5, 0.5]) == 0

        @test is_ccw(triang, [0,0,1])

        square = [0 1 1 0; 0 0 1 1; 0 0 0 0]
        @test earcut3d(square) == [[3,4,1], [1,2,3]]
        @test inpolygon3d(square, center_of_mass(square)) == 1
        @test inpolygon3d(square, [1,1,1]) == 0
        @test inpolygon3d(square, [1/2, 0,0]) == -1
        @test inpolygon3d(square, [0,0,0]) == -1
    end

    @testset "Polyhedron.jl" begin
        # Polyhedron construction, TODO: get_... and set_... testen
        @test_throws "dimension 2" Polyhedron(rand(Float64, (3,4)), [[1,2], [2,3], [3,4], [4,1]], [[1,2,3,4]])
        
        tetrahedron1 = Polyhedron([0 1 0 0; 0 0 1 0; 0 0 0 1], [[1,2], [2,3], [3,1], [4,2], [4,3], [4,1]], [[1,2,3], [2,3,4], [1,2,4], [3,4,1]])
        tetrahedron2 = Polyhedron([0 1 0 0; 0 0 1 0; 1 1 1 2], [[1,2], [2,3], [3,1], [4,2], [4,3], [4,1]], [[1,2,3], [2,3,4], [1,2,4], [3,4,1]])
        @test iscongruent(tetrahedron1, tetrahedron2)

        butterfly = Polyhedron(rand(Float64, (3,4)), [[1,2], [2,3], [3,1], [2,4], [3,4]], [[1,2,3], [2,3,4]])
        @test isadjacent(butterfly, [2,3], [1,2,3])
        @test isadjacent(butterfly, [1,2,3], [2,3,4])
        @test Set(Set.(adjfacets(butterfly, [2,3]))) == Set(Set.([[1,2,3], [2,3,4]]))
        @test isincident(1, [1,2,3])
        @test isincident(1, [1,2])
        @test incfacets(butterfly, [1,2,4]) == []
        @test Set.(incfacets(butterfly, [1,2])) == Set.([[1,2,3]])
        @test Set.(incfacets(butterfly, 1)) == Set.([[1,2,3]])
        @test Set(Set.(incedges(butterfly, 1))) == Set(Set.([[1,2], [3,1]]))
        @test incedges(butterfly, [1,2,3]) == [[1,2], [2,3], [3,1]]
        @test inpolyhedron(center_of_mass(get_verts(butterfly)[:, 1:3]), butterfly) == -1
        
        tetrahedron = Polyhedron([0 2 0 0; 0 0 2 0; 0 0 0 2], [[1,2], [2,3], [3,1], [4,2], [4,3], [4,1]], [[1,2,3], [4,3,2], [4,2,1], [3,4,1]])
        @test inpolyhedron(center_of_mass(get_verts(tetrahedron)), tetrahedron) == 1
        @test inpolyhedron([2,2,2], tetrahedron) == 0

        @test vol_signed(tetrahedron) < 0
        set_facets!(tetrahedron, reverse.(get_facets(tetrahedron)))
        @test vol_signed(tetrahedron) > 0

        for i in 1:10
            coords = hcat([0 1 1 0; 0 0 1 1; 0 0 0 0], rand(Float64, 3))
            h = coords[3, 5]
            pyramide = Polyhedron(coords, [[1,2], [2,3], [3,4], [4,1], [1,5], [2,5], [3,5], [4,5]], [[1,2,3,4], [1,2,5], [2,3,5], [3,4,5], [1,4,5]])
            @test vol(pyramide) ≈ 1/3 * h

            d = rand(Float64)
            coords = d * [0 1 1 0; 0 0 1 1] + rand(Float64) * ones((2,4))
            square = Polyhedron(coords, [[1,2], [2,3], [3,4], [4,1]], [[1,2,3,4]])
            @test vol(square) ≈ d^2
        end
    end

    @testset "merging.jl" begin
        function nprism(n::Integer)
            alpha = 1/dist([1,0], [cos(2*pi/n), sin(2*pi/n)])
            verts = vcat([[alpha * cos(2*pi*k/n), alpha * sin(2*pi*k/n), 0] for k in 0:n-1], [[alpha * cos(2*pi*k/n), alpha * sin(2*pi*k/n), 1] for k in 0:n-1])
            facets = vcat([[1:n...]], [[(n+1):(2*n)...]], [[k, k+1, k+n+1, k+n] for k in 1:(n-1)], [[n,1,n+1,2*n]])
            edges = vcat([[k, k+1] for k in 1:(n-1)], [[n,1]], [[n+k, n+k+1] for k in 1:(n-1)], [[2*n, n+1]], [[k, k+n] for k in 1:n])
            return Polyhedron(verts, edges, facets)
        end

        cube1 = nprism(4)
        cube2 = nprism(4)

        @test_throws DimensionMismatch merge(cube1, cube2, [[1,2,3,4]], [[1,2,3,4], [5,6,2,1]])
        @test_throws "consist of facets" merge(cube1, cube2, [[1,2,3,4]], [[1,2,3,5]])
        @test_throws "well defined" merge(cube1, cube2, [[1,2,3,4], [1,2,3,4]], [[1,2,3,4], [5,6,2,1]])
        poly = merge(cube1, cube2, [[1,2,3,4]], [[1,2,3,4]])
        @test size(get_verts(poly))[2] == size(get_verts(cube1))[2] + size(get_verts(cube2))[2] - 4
        @test sort(unique(vcat(get_edges(poly)...))) == collect(1:12)
        @test length(get_edges(poly)) == length(get_edges(cube1)) + length(get_edges(cube2)) - 4
        @test sort(unique(vcat(get_facets(poly)...))) == collect(1:12)
        @test length(get_facets(poly)) == length(get_facets(cube1)) + length(get_facets(cube2)) - 2
        @test_throws "cannot be merged" merge(cube1, cube2, [[1,2,3,4]], [[2,1,3,4]])

        set_verts!(cube2, 3 * get_verts(cube2))
        @test_throws "cannot be merged" merge(cube1, cube2, [[1,2,3,4]], [[1,2,3,4]])
    end

    @testset "decomposition.jl" begin
        function nprism(n::Integer)
            alpha = 1/dist([1,0], [cos(2*pi/n), sin(2*pi/n)])
            verts = vcat([[alpha * cos(2*pi*k/n), alpha * sin(2*pi*k/n), 0] for k in 0:n-1], [[alpha * cos(2*pi*k/n), alpha * sin(2*pi*k/n), 1] for k in 0:n-1])
            facets = vcat([[1:n...]], [[(n+1):(2*n)...]], [[k, k+1, k+n+1, k+n] for k in 1:(n-1)], [[n,1,n+1,2*n]])
            edges = vcat([[k, k+1] for k in 1:(n-1)], [[n,1]], [[n+k, n+k+1] for k in 1:(n-1)], [[2*n, n+1]], [[k, k+n] for k in 1:n])
            return Polyhedron(verts, edges, facets)
        end

        for n in 3:10
            @test length(get_facets(triangulate(nprism(n)))) == 2*n + 2*(n-2)
        end

        tetrahedron = Polyhedron([0 2 0 0; 0 0 2 0; 0 0 0 2], [[1,2], [2,3], [3,1], [4,2], [4,3], [4,1]], [[1,2,3], [4,3,2], [4,2,1], [3,4,1]])
        @test outward_normal(tetrahedron, [1,2,3]) == [0,0,-1]
        set_facets!(tetrahedron, reverse.(get_facets(tetrahedron)))
        @test outward_normal(tetrahedron, [3,2,1]) == [0,0,-1]
        
        @test isconvex(tetrahedron) == true

        function tiblock(n1::Integer, n2::Integer, nmerges::Integer; atol::Real = 1e-8)
            @assert mod(n1, nmerges) == 0 "Number of blocks on the side needs to divide n1."
            
            if n2 == 3
                sideelem = nprism(3)
                merge!(sideelem, nprism(3), [[3,1,6,4]], [[1,2,4,5]], atol = atol)
                merge!(sideelem, nprism(3), [[2,3,5,6]], [[1,2,4,5]], atol = atol)
            elseif n2 == 4
                sideelem = nprism(4)
                merge!(sideelem, nprism(4), [[2,3,6,7]], [[1,2,5,6]], atol = atol)
                merge!(sideelem, nprism(4), [[4,1,8,5]], [[1,2,5,6]], atol = atol)
            else
                sideelem = nprism(n2)
            end
            
            block = nprism(n1)
            for k in 3:Int(n1 / nmerges):length(get_facets(nprism(n1)))
                merge!(block, sideelem, [get_facets(nprism(n1))[k]], [[n2+1,1,2, n2+2]], atol = atol)
            end
            
            return block
        end
        poly = tiblock(6,3,2)
        edgetypes = [edgetype(poly, e) for e in get_edges(poly)]
        @test count(e -> e == "flat", edgetypes) == 10
        @test count(e -> e == "convex", edgetypes) == 24
        @test count(e -> e == "concave", edgetypes) == 14
        @test isconvex(poly) == false
        @test length(get_facets(flattenfacets(poly))) == 16
        @test length(get_edges(flattenfacets(poly))) == 38
    end

    @testset "interlocking.jl" begin
        
    end

    @testset "SimplicialSurface.jl" begin
        vertices=1.0*[[1,0,1],[ 1, 1, 0 ], [ 2, 1, 1 ],[ 1, 1, 2 ], [ 0, 1, 1 ],[ 1, 2, 1 ]]
        facets=[ [ 1, 2, 3 ], [ 2, 5, 6 ], [ 1, 2, 5 ], [ 2, 3, 6 ], [ 1, 4, 5 ], 
        [ 3, 4, 6 ], [ 1, 3, 4 ], [ 4, 5, 6 ] ]
        edges=[ [ 1, 2 ], [ 1, 3 ], [ 1, 4 ], [ 1, 5 ], [ 2, 3 ], [ 2, 5 ], [ 2, 6 ], 
        [ 3, 4 ], [ 3, 6 ], [ 4, 5 ], [ 4, 6 ], [ 5, 6 ] ]
        Octahedron=Polyhedron(vertices,edges,facets)

        surf = SimplicialSurface(Octahedron)
        @test length(boundary(surf)) == 0
        @test length(get_verts(surf)) == length(get_verts(edgeturn(surf, [1,2])))
        @test length(get_edges(surf)) == length(get_edges(edgeturn(surf, [1,2])))
        @test length(get_facets(surf)) == length(get_facets(edgeturn(surf, [1,2])))
        @test length(Base.intersect([[1,2], [2,1]], get_edges(edgeturn(surf, [1,2])))) == 0

        Tetrahedron=Polyhedron([[0.5,0,0],[-0.5,0.,0.],[0.,0.5,sqrt(3.)/2.],[0.,-0.5,sqrt(3.)/2.]],
            [ [ 1, 2 ], [1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ], [ 3, 4 ] ],
            [ [ 1, 2, 3 ], [ 1, 2, 4 ], [ 2, 3, 4 ], [ 1, 3, 4 ] ])
        surf = SimplicialSurface(Tetrahedron)
        @test_throws "not turnable" edgeturn(surf, [1,2])

        set_facets!(surf, [[1,2,3], [1,2,4], [2,3,4]])
        @test length(boundary(surf)) == 3
    end
end

test()