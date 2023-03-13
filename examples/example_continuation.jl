using Plots
using NS2DBenchmarkSolver

nbcells = (32, 32)
reynolds = 100 # already after pitch fork bifurcation starts
reynoldsrange = reynolds:reynolds+200

mesh = SpectralMesh2D(nbcells)
probl = Cavity4Sided(mesh, reynolds)

k0 = 10
bcfpos(x) = @. ((exp(k0 * (x - 1)) - 1) * (exp(-k0 * (x + 1)) - 1))^2
bcfneg(x) = -bcfpos(x)

setBC2D(probl, bcfneg, bcfpos, bcfpos, bcfneg)

# Create initiual guess for a stable solution
Ψinit = NS2DBenchmarkSolver.calculateinitialguess(probl, Δt=1, nbtimesteps=130)
println("Initial guess calculated")

nx = probl.mesh.xnbcells
ny = probl.mesh.ynbcells
nxc = Int(floor(probl.mesh.xnbcells/2)+1)
nyc = Int(floor(probl.mesh.ynbcells/2)+1)

ψcenter = Vector{Float64}()

for re in reynoldsrange
    probl.reynolds = re

    sol = NS2DBenchmarkSolver.solve(probl, Ψinit)

    println("Reynolds $re : $(sol.iter) iterations")

    center = [Ψinit[nxc, nyc]]
    global ψcenter = append!(ψcenter, center)

    global Ψinit = sol.vals
end

plot(reynoldsrange,
    ψcenter,
    title="Bifurcation",
    xlabel="Re",
    ylabel="Ψcenter",
    label=false
)
