using FourSidedCavityFlow
using NonlinearSolve, DelimitedFiles, Plots
const CF = FourSidedCavityFlow

N = 32
Re = 100
p = CF.setup_struct(N, Re)

# initial condition
Ψ0 = zeros(N + 1, N + 1)
u0 = Ψ0[3:(N - 1), 3:(N - 1)][:]

# solver
probN = NonlinearProblem(CF.f!, u0, p)
solver = @time solve(probN, NewtonRaphson(autodiff=false), reltol = 1e-11)

Ψsol = CF.constructBC(solver.u, p)
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
