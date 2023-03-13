using Plots
using NS2DBenchmarkSolver

# Function implementing timestepping and creating a gif
# of the time evolution of the streamfunction for the cavity problem
function timestep(probl::Cavity4Sided, Ψinit, Δt, nbtimesteps)
    plt = heatmap(
        xlim=(-1,1),
        ylim=(-1,1),
        title="Streamfunction Ψ on $nx x $ny",
        titlefontsize=10,
        xtickfontsize=8,
        aspect_ratio=1,
        legend=false,
    )
    plt2 = scatter(
        xlim=(0,Δt*nbtimesteps),
        ylim=(-0.25,0.25),
        title="Ψcenter value",
        titlefontsize=10,
        xtickfontsize=8,
        ytickfontsize=8,
        legend=false,
        markercolor=:blue,
    )

    Ψold = Ψinit
    ψcenter = [Ψold[nxc, nyc]]
    time = 0

    anim = @animate for step in 1:nbtimesteps
        println(step)
        plot(heatmap!(plt, probl.mesh.xnodes, probl.mesh.ynodes, Ψold),
            scatter!(plt2, time, ψcenter, markercolor = :blue),
            layout=grid(1,2, widths=(0.5,0.5)), size=(600,300)
        )

        ψint = vec(Ψold[3:nx-1, 3:ny-1])

        function ftimestep(ψint) 
            return NS2DBenchmarkSolver.rhstime(probl, Δt, Ψold, ψint)
        end

        ψintold = vec(Ψold[3:nx-1, 3:ny-1])
        ψint, _, _, _ = NS2DBenchmarkSolver.newton(ftimestep, ψintold)

        Ψint = reshape(ψint, (nx-3,nx-3))
        Ψold = NS2DBenchmarkSolver.constructΨboundary(Ψint, probl.mesh.diffx1mat, probl.mesh.diffy1mat, probl.bcxmin, probl.bcxmax, probl.bcymin, probl.bcymax)

        ψcenter = [Ψold[nxc, nyc]]
        time = [(Δt * step)]
    end

    gif(anim, "example_timestepping.gif", fps=15)
end


# Example problem
nbcells = (32, 32)
reynolds = 100

mesh = SpectralMesh2D(nbcells)
probl = Cavity4Sided(mesh, reynolds)

k0 = 10
bcfpos(x) = @. ((exp(k0 * (x - 1)) - 1) * (exp(-k0 * (x + 1)) - 1))^2
bcfneg(x) = -bcfpos(x)

setBC2D(probl, bcfneg, bcfpos, bcfpos, bcfneg)

Δt = 1
nbtimesteps = 140

nx = probl.mesh.xnbcells
ny = probl.mesh.ynbcells
nxc = Int(floor(probl.mesh.xnbcells/2)+1)
nyc = Int(floor(probl.mesh.ynbcells/2)+1)

# Convergence to unstable solution
# Ψinit = zeros((nx+1),(ny+1))

# Convergence to one of the two possible assymetric solutions
Ψinit = 1e-3*randn((nx+1),(ny+1))

timestep(probl, Ψinit, Δt, nbtimesteps)
