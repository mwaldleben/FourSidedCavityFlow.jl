using FourSidedCavityFlow
using BifurcationKit, DelimitedFiles, Plots
const CF = FourSidedCavityFlow
const BK = BifurcationKit

N = 64
Re = 100
p = CF.setup_struct(N, Re)

# initial condition
Ψ0 = zeros(N + 1, N + 1)
u0 = Ψ0[3:(N - 1), 3:(N - 1)][:]

# solver
prob = BK.BifurcationProblem(CF.f, u0, p, (@lens _.Re))
opt_newton = BK.NewtonPar(tol = 1e-8, verbose = true, maxIter = 200)
sol = @time BK.newton(prob, opt_newton)

Ψsol = CF.constructBC(sol.u, p)
display(Ψsol)

plt = contourf(reverse(p.params.nodes),
               reverse(p.params.nodes),
               Ψsol',
               xlim = (-1, 1),
               ylim = (-1, 1),
               aspect_ratio = 1,
               axis = ([], false),
               legend = false,
               c = :davos)
savefig(plt, "psi.png")
