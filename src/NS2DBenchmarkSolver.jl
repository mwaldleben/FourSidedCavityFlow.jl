module NS2DBenchmarkSolver

using LinearAlgebra

include("typedefs.jl")
include("differentiation.jl")
include("solve.jl")
include("boundarycondition.jl")
include("spectralmesh.jl")
include("spectralproblem.jl")

export Solution1D
export SpectralMesh1D
export SpectralMesh2D
export Example1D
export Cavity4Sided
export setNeumannBC2D
export calculateΨinitial
export solve

end
