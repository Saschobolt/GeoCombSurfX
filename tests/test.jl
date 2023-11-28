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
        tetrahedron = Polyhedron([0 2 0 0; 0 0 2 0; 0 0 0 2], [[1,2], [2,3], [3,1], [4,2], [4,3], [4,1]], [[1,2,3], [2,3,4], [1,2,4], [3,4,1]])
        @test inpolyhedron(center_of_mass(get_verts(tetrahedron)), tetrahedron) == 1
        @test inpolyhedron([2,2,2], tetrahedron) == 0
    end
end

test()