module FourSidedCavityFlow

using LinearAlgebra
using FiniteDiff
using UnPack
using Printf

include("differentiation.jl")
include("newtonraphson.jl")
include("foursidedcavity.jl")
include("steadystate.jl")
include("timestepping.jl")
include("continuation.jl")
include("linearstability.jl")

end
