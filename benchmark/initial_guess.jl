import CavityFlow as CF
using DelimitedFiles
using Random
gr()

Random.seed!(1234)

# Timestepping to get an intial guess for the benchmark
n = 64
Re = 100

p = CF.setup_params(n, Re)

Δt = 1
steps = 140
Ψ0 = 1e-3 * randn((n + 1), (n + 1))

Ψ = CF.solve_timestepping(Ψ0, p, Δt, steps)

writedlm("$(n)x$(n)_initial_guess_Re$(Re).txt", Ψ)
