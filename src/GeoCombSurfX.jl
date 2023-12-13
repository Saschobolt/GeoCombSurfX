include("affine_geometry.jl")
include("Framework.jl")
include("Polyhedron.jl")
include("merging.jl")
include("decomposition.jl")
include("read.jl")
include("plotting.jl")
# include("examples.jl")

# ##################################################################################
# ################# Module: GeoCombSurfX
# ##################################################################################
# module GeoCombSurfX
# include("affine_geometry.jl")
# export indcols_indices, indcols, colspace, affinebasis_indices, affinebasis, affinedim, affinemap, rigidmap, normalvec, sqdist, signedangle3d_right, signedangle3d_left, center_of_mass

# include("Framework.jl")
# export AbstractSimpleGraph, Graph, get_verts, get_edges, set_verts!, set_edges!, no_concomponents, AbstractEmbeddedGraph, Framework, dimension, display

# include("Polyhedron.jl")
# export AbstractPolyhedron, Polyhedron, get_facets, set_facets!, is_congruent, display, isadjacent, adjfacets, isincident, incfacets, incedges, inpolyhedron

# include("merging.jl")
# export merge!, merge

# include("decomposition.jl")
# export triangulate!, triangulate, outward_normal, isflatedge, edgetype, flattenfacets!, flattenfacets, isturnable, isconvex

# include("read.jl")
# export readPolyhedronFromSTL, readAssemblyFromSTL, readPolyhedronFromTextfile, reArrangePolyhedron

# # include("plotting.jl")
# # export trace_polyhedron, plot

# include("examples.jl")
# export nprism, Octahedron, Cube, Tetrahedron, Icosahedron, Dodecahedron

# ##################################################################################
# #################  Combinatorics
# ##################################################################################
# # include("Combinatorics//combinatorics.jl")

# ##################################################################################
# #################  Interlocking
# ##################################################################################
# # include("examples.jl")
# export tiblock, assembly1, assembly2, assembly3, assembly4

# include("Interlocking/interlocking.jl")
# export contacts, titest

# ##################################################################################
# #################  Rigidity
# ##################################################################################
# # include("examples.jl")
# export bricard_octahedron

# ##################################################################################
# #################  Embeding
# ##################################################################################

# end;