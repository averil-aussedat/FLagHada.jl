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
    include(theexample * "/structure.jl")
    include(theexample * "/geometry.jl")
    include(theexample * "/mesh.jl")
    include(theexample * "/interface.jl")
    include(theexample * "/data.jl")
end

end
