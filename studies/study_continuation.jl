using FourSidedCavityFlow
using Plots
const CF = FourSidedCavityFlow

Re_start = 100 
p.params.Re = Re_start

println("--- Continuation: $(n)x$(n) ---")
println("Calculate unstable solution at Re_start = $Re_start :")
steps = 140
Δt = 1
Ψi = 1e-3 * randn((n + 1), (n + 1))
Ψt = @time CF.timestepping(Ψi, p, Δt, steps)
Ψs, iter, tol = CF.steadystate(Ψt, p)

ΔRe = 1
steps = 840

println("Arclength continuation, ΔRe = $ΔRe, steps = $steps :")
@time CF.continuation_arclength(Ψs, p, Re_start, ΔRe, steps)
