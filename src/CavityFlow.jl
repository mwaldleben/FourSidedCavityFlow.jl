module CavityFlow

using LinearAlgebra

include("typedefs.jl")
include("differentiation.jl")
include("solve.jl")
include("boundarycondition.jl")
include("rectangularmesh.jl")
include("cavityproblem.jl")

export ChebyshevMesh
export Cavity4Sided
export setNeumannBC
export calculateΨinitial
export solve

end
