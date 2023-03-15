using Plots
using NS2DBenchmarkSolver

n = (32, 32)
reynoldsrange = 100:-1:10 # we want to go backwards

mesh = SpectralMesh2D(n)
probl = Cavity4Sided(mesh, reynolds)

k0 = 10
bcfunc(x) = @. ((exp(k0*(x-1)) - 1) * (exp(-k0*(x+1)) - 1))^2
bcfuncneg(x) = -bcfunc(x)

setNeumannBC2D(probl, bcfuncneg, bcfunc, bcfunc, bcfuncneg)

# Create initiual guess for a stable solution
println("Calculating initial guess...")
Ψinitial = NS2DBenchmarkSolver.calculateΨinitial(probl, Δt=1, nbtimesteps=140)
contourf(reverse(mesh.xnodes), reverse(mesh.ynodes), Ψinitial')

nx = probl.mesh.nx
ny = probl.mesh.ny
nxc = Int(floor(probl.mesh.nx/2)+1)
nyc = Int(floor(probl.mesh.ny/2)+1)

ψcenter = Vector{Float64}()
for re in reynoldsrange
    probl.reynolds = re

    global Ψinitial, iter, _, _ = solve(probl, Ψinitial)

    println("Reynolds $re : $(iter) iterations")

    center = [Ψinitial[nxc, nyc]]
    global ψcenter = append!(ψcenter, center)
end

plot(reverse(reynoldsrange),
    reverse(ψcenter),
    title="Bifurcation",
    xlabel="Re",
    ylabel="Ψcenter",
    label=false
)
