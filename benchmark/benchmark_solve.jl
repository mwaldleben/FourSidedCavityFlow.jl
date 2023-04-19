import CavityFlow as CF
using BenchmarkTools
using DelimitedFiles
using FiniteDiff
using Random

Random.seed!(1234)

# Example of the cavity flow to find a steady state solution 
n = 32
Re = 100

p = CF.setup_params(n, Re)

filename =  "$(n)x$(n)_initial_guess_Re$(Re).txt"
filepath = joinpath(@__DIR__, filename)

Ψ0 = readdlm(filepath)

u0 = reshape(Ψ0[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
fu0 = similar(u0)

##### Benchmark function evaluation #####
display(@benchmark CF.f!(fu0, u0, p))

##### Benchmark Jacobian #####
dim = size(u0, 1)
J = zeros(dim, dim)
cache = FiniteDiff.JacobianCache(u0)
f1!(fu0, u0) = CF.f!(fu0, u0, p)

display(@benchmark CF.jacobian!(J, f1!, u0))
display(@benchmark FiniteDiff.finite_difference_jacobian!(J, f1!, u0, cache))

##### Benchmark solve #####
display(@benchmark CF.solve_steadystate(Ψ0, p))
