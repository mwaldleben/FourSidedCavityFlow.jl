module FourSidedCavityFlow

using LinearAlgebra
using FiniteDiff
using UnPack

include("differentiation.jl")
include("newtonraphson.jl")
include("foursidedcavity.jl")
include("steadystate.jl")
include("timestepping.jl")
include("linearstability.jl")

end
