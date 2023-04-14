import CavityFlow as CF
using Plots
gr()

# Example of the cavity flow to find a steady state solution 
n = 64
Re = 30

p = CF.setup_params(n, Re)
Ψ0 = zeros(n + 1, n + 1)

Ψ, iter, tol = @time CF.solve_steadystate(Ψ0, p)
println("Iterations: $iter")

# Reverse solution as non ascending intervals are not supported by Plots
contourf(reverse(p.nodes), reverse(p.nodes), Ψ', aspect_ratio = 1)
