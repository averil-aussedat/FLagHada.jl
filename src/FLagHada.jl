module FLagHada

using LinearAlgebra
using Parameters
using DocStringExtensions 
using ProgressBars
using Polynomials
using Plots
using Colors
using LazySets, Polyhedra

include("structure.jl")
include("scheme.jl")
include("myutils.jl")
include("interface.jl")

for theexample in ["ellipses", "euclideans", "networks", "gluings"]
    include(theexample * "/structure.jl") # must provide Point, Domain, Mesh
    include(theexample * "/geometry.jl") # must provide get_distance, get_convex, get_exponential
    include(theexample * "/mesh.jl") # must provide get_mesh, search_in_mesh, get_mesh_convex 
    include(theexample * "/interface.jl") # plotting functions, no requirement
    include(theexample * "/data.jl") # constructor (usually get_[type], sometimes a plotter), get_Scalar (cost), get_Dynamic, get_Analytical
end

end
