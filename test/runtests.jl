import FourSidedCavityFlow
const CF = FourSidedCavityFlow

using LinearAlgebra
using FileIO
using Test

include("test_differentiation.jl")
include("test_newtonraphson.jl")
include("test_foursidedcavity.jl")
include("test_steadystate.jl")
include("test_timestepping.jl")
include("test_continuation.jl")
include("test_linearstability.jl")
