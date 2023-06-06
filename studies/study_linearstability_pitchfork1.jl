import FourSidedCavityFlow as CF
using DelimitedFiles

# Read a steady state Ψ generated with "study_steadystate.jl"
n = 32
Re0 = 66
p = CF.setup_params(n, Re0)
name = "sym"

title = "--- Linear stability analysis: $(n)x$(n) ---"
println(title)

filename = "results/$(n)x$(n)_$(name)_Re$(Re0).txt"
filepath = joinpath(@__DIR__, filename)
println("... reading Ψstart from $(filename)")
Ψstart = readdlm(filepath)

u0 = Ψstart[3:(n - 1), 3:(n - 1)][:]

println("1D Newton:")
Re, iter, tol = CF.newton1D_for_linearstability(Re0, u0, p; tolmax = 1e-12,
                                                maxiter = 20, verbose = true)

println("Re: $Re")
println(repeat('-', length(title)))

