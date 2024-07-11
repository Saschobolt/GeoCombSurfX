# include("affine_geometry.jl")
# include("Framework.jl")
# include("Polyhedron.jl")
# include("merging.jl")
# include("decomposition.jl")
# # include("read.jl")
# include("plotting.jl")
# # include("examples.jl")

module GeoCombSurfX

using LinearAlgebra
using StaticArrays
using PlotlyJS
using Colors
using JuMP

import HiGHS
import Polyhedra
import Graphs
import Graphs.SimpleGraphs, Graphs.connected_components
import Base.Multimedia.display
import PlotlyJS.plot

export Ray, Plane
export indcols_indices, indcols, colspace, affinebasis_indices, affinebasis, affinedim, affinemap, rigidmap, normalvec, dist, sqdist, signedangle2d, signedangle3d_right, signedangle3d_left, center_of_mass, bounding_sphere

export to_xyplane_map, to_xyplane, intriang, is_ccw, earcut3d, inpolygon3d

export AbatractEmbeddedGraph
export Framework
export get_verts, get_edges, adjacency_matrix, get_edges, get_verts, set_verts!, dimension, edge_lengths, henneberg_extension!, henneberg_extension

export AbstractPolyhedron, AbstractCombPolyhedron, AbstractEmbOrCombPolyhedron
export HalfEdge, Polyhedron
export halfedge, halfedges, set_halfedges!, get_verts, set_verts!, set_edges!, get_facets, set_facets!, iscongruent, dimension, orient_facets!, orient_facets, isadjacent, adjfacets, isincident, incfacets, incedges, facet_adjacency, boundary, removefacet!, removefacet, removeedge!, removeedge, vol_signed, orient_facets_ccw!, orient_facets_ccw, vol

export AbstractCombSimplicialSurface, AbstractSimplicialSurface, AbstractEmbOrCombSimplicialSurface, AbstractColoredSimplicialSurface
export CombSimplicialSurface, SimplicialSurface, ColoredSimplicialSurface
export vertex_degree, characterictic, iscactus, is_SimplicialSurface, set_facets!, remove_vertex!, insert_butterfly!, insert_butterfly, random_simplsphere, random_emb_simplsphere, edgeturn!, edgeturn, remove_tetrahedron!, append_tetrahedron!, random_cactus, random_emb_cactus
export colors, color, colortype, congruencetypes, is_congcolored, edge_type, is_tamecolored

export merge!, merge # merge and merge! are also exported by Base

export triangulate!, triangulate, outward_normal, isflatedge, edgetype, remove_flatedge!, remove_edge, flattenfacets!, flattenfacets, isconvex

export plot

export titest

export rigidity_matrix, basis_inf_motions, is_infrigid, basis_inf_flex, index, is_genrigid

include("affine_geometry.jl")
include("polygonal_geometry.jl")

abstract type AbstractEmbeddedGraph{S<:Real,T<:Integer} end

include("Framework.jl")

abstract type AbstractPolyhedron{S<:Real,T<:Integer} <: AbstractEmbeddedGraph{S,T} end
abstract type AbstractCombPolyhedron{T<:Integer} end
AbstractEmbOrCombPolyhedron = Union{AbstractPolyhedron,AbstractCombPolyhedron}

include("Polyhedron.jl")

abstract type AbstractCombSimplicialSurface{T<:Integer} <: AbstractCombPolyhedron{T} end
abstract type AbstractSimplicialSurface{S<:Real,T<:Integer} <: AbstractPolyhedron{S,T} end
AbstractEmbOrCombSimplicialSurface{S<:Real,T<:Integer} = Union{AbstractSimplicialSurface{S,T},AbstractCombSimplicialSurface{T}}
abstract type AbstractColoredSimplicialSurface{T<:Integer} <: AbstractCombSimplicialSurface{T} end


include("SimplicialSurface.jl")
include("merging.jl")
include("decomposition.jl")
include("plotting.jl")

include("Interlocking/interlocking.jl")
include("Rigidity/infinitesimal_rigidity.jl")


end # module
