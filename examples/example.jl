using Plots
include("../src/NS2DBenchmarkSolver.jl")

nbcells = (32, 32)

mesh = NS2DBenchmarkSolver.SpectralMesh2D(nbcells)

reynolds = 50
probl = NS2DBenchmarkSolver.Cavity4Sided(mesh, reynolds)

k0 = 10
bcfpos(x) = @. ((exp(k0 * (x - 1)) - 1) * (exp(-k0 * (x + 1)) - 1))^2
bcfneg(x) = -bcfpos(x)

NS2DBenchmarkSolver.setBC2D(probl, bcfneg, bcfpos, bcfpos, bcfneg)

sol = NS2DBenchmarkSolver.solve(probl)

NS2DBenchmarkSolver.print(sol)
plot(sol)
