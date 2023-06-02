module FourSidedCavityFlow

using LinearAlgebra
using FiniteDiff
using PreallocationTools
using UnPack
using DelimitedFiles
using Printf
using UnPack

include("differentiation.jl")
include("newtonraphson.jl")
include("foursidedcavity.jl")
include("steadystate.jl")
include("timestepping.jl")
include("continuation.jl")
include("linearstability.jl")

end
