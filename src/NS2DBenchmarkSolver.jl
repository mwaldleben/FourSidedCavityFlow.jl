module NS2DBenchmarkSolver

using LinearAlgebra

include("typedefs.jl")
include("differentiation.jl")
include("boundarycondition.jl")
include("solve.jl")
include("spectralmesh.jl")
include("spectralproblem.jl")
include("solution.jl")

export BCDirichlet1D 
export BCNeumann1D
export BCDirichlet2D 
export BCNeumann2D
export Solution1D
export SpectralMesh1D
export SpectralMesh2D
export ExamplePoisson1D
export Cavity4Sided
export setBC2D
export solve

end
