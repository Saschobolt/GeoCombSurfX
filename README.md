# GeoCombSurfX

This package aims to facilitate the exploration and manipulation of piecewise linear surfaces, particularly simplicial surfaces.

## Basic Usage
A simplicial surface is a piecewise linear surface with triangular faces. To construct a simplicial surface using this package, provide the combinatorial structure:
```Julia
# construct tetrahedron as a cominatorial structure.
tetra = CombSimplicialSurface(verts = [1,2,3,4], edges = [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], facets = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]])
```
The constructed combinatorial simplicial surface is a purely combinatorial object. To construct an embedded tetrahedron, construct it as a SimplicialSurface and provide the vertex coordinates:
```Julia
# construct tetrahedron with embedded vertices.
tetra = SimplicialSurface(verts = [0 1 0 0; 0 0 1 0; 0 0 0 1], edges = [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], facets = [[1,2,3], [1,2,4], [1,3,4], [2,3,4]])

# orient the facets of the tetrahedron counter clockwise wrt the outward facing normal vectors.
orient_facets_ccw!(tetra)
```