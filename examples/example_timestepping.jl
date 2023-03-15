using Plots
using CavityFlow

# Function implementing timestepping and creating a gif
# of the time evolution of the streamfunction for the cavity problem
function timestep(probl::Cavity4Sided, Ψinit, Δt, nbtimesteps)
    plt = contourf(
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
        plot(contourf!(plt, reverse(probl.mesh.xnodes), reverse(probl.mesh.ynodes), Ψold'),
            scatter!(plt2, time, ψcenter, markercolor = :blue),
            layout=grid(1,2, widths=(0.5,0.5)), size=(600,300)
        )

        ψi = vec(Ψold[3:nx-1, 3:ny-1])

        function ftimestep(ψint) 
            return CavityFlow.rhstime(probl, Δt, Ψold, ψint)
        end

        ψiold = vec(Ψold[3:nx-1, 3:ny-1])
        ψi, _, _, _ = CavityFlow.newton(ftimestep, ψiold)

        Ψi = reshape(ψi, (nx-3,nx-3))
        Ψold = CavityFlow.constructΨ(probl, Ψi)

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
bcfunc(x) = @. ((exp(k0*(x-1)) - 1) * (exp(-k0*(x+1)) - 1))^2
bcfuncneg(x) = -bcfunc(x)

setNeumannBC2D(probl, bcfuncneg, bcfunc, bcfunc, bcfuncneg)

Δt = 1
nbtimesteps = 140

nx = probl.mesh.nx
ny = probl.mesh.ny
nxc = Int(floor(nx/2)+1)
nyc = Int(floor(ny/2)+1)

# Convergence to unstable solution
# Ψinit = zeros((nx+1),(ny+1))

# Convergence to one of the two possible assymetric solutions
Ψinitial = 1e-3*randn((nx+1),(ny+1))

timestep(probl, Ψinitial, Δt, nbtimesteps)
