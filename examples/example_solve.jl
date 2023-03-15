using Plots
using CavityFlow

# Example of the cavity flow to find a steady state solution 
# using an initial guess generated from time evolution
n = (32, 32)
reynolds = 100

mesh = SpectralMesh2D(n)
probl = Cavity4Sided(mesh, reynolds)

k0 = 10
bcfunc(x) = @. ((exp(k0*(x-1)) - 1) * (exp(-k0*(x+1)) - 1))^2
bcfuncneg(x) = -bcfunc(x)
setNeumannBC2D(probl, bcfuncneg, bcfunc, bcfunc, bcfuncneg)

Ψinitial = calculateΨinitial(probl; nbtimesteps=1)
Ψ, iter, tol, isconverged = solve(probl, Ψinitial)
println("Converged: $isconverged, iterations: $iter")

# Reverse solution as non ascending intervals are not supported by Plots
contourf(reverse(mesh.xnodes), reverse(mesh.ynodes), Ψ')
