using FourSidedCavityFlow
using Random
const CF = FourSidedCavityFlow

Random.seed!(1234)

n = 32
Re = 1 # dummy value
p = CF.setup_struct(n, 1)

# include("study_continuation.jl")
# include("study_converge_psis.jl")
include("study_converge_psis_bifurcation.jl")
include("study_bifurcation_diag.jl")
