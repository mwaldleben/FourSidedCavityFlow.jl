module FourSidedCavityFlow

using LinearAlgebra
using FiniteDiff
using DelimitedFiles
using FileIO
using JLD2
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
