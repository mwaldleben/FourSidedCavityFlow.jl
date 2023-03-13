using Plots
using NS2DBenchmarkSolver

# Example of the cavity flow to find a steady state solution 
# using an initial guess generated from time evolution
nbcells = (32, 32)
reynolds = 100

mesh = NS2DBenchmarkSolver.SpectralMesh2D(nbcells)
probl = NS2DBenchmarkSolver.Cavity4Sided(mesh, reynolds)

k0 = 10
bcfpos(x) = @. ((exp(k0 * (x - 1)) - 1) * (exp(-k0 * (x + 1)) - 1))^2
bcfneg(x) = -bcfpos(x)
NS2DBenchmarkSolver.setBC2D(probl, bcfneg, bcfpos, bcfpos, bcfneg)

Ψinit = NS2DBenchmarkSolver.calculateinitialguess(probl; nbtimesteps=150)
sol = NS2DBenchmarkSolver.solve(probl, Ψinit)

NS2DBenchmarkSolver.print(sol)
plot(sol)
