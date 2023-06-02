using FourSidedCavityFlow
using Plots
using Random
const CF = FourSidedCavityFlow

Random.seed!(1234)

n = 32
Re_start = 100
p = CF.setup_struct(n, Re_start)

println("--- Continuation: $(n)x$(n), Re_start = $Re_start ---")

println("Calculate unstable solution at Re_start = $Re_start...")
steps = 140
Δt = 1
Ψi = 1e-3 * randn((n + 1), (n + 1))
Ψt = @time CF.timestepping(Ψi, p, Δt, steps; verbose=true)
Ψs, iter, tol = CF.steadystate(Ψt, p)

println("Arclength continuation...")
ΔRe = 1
steps = 830
CF.continuation_arclength(Ψs, p, Re_start, ΔRe, steps)
