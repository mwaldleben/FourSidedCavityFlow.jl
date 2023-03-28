using Plots; gr()
using CavityFlow
const CF = CavityFlow

# Timestepping and creating a gif of the time evolution of the streamfunction 
# for the cavity problem
n = 32 
Re = 100

p = CF.setup_params(n, Re)

Δt = 1
steps = 140

nxc = Int(floor(n/2)+1)
nyc = Int(floor(n/2)+1)

# Convergence to one of the two possible assymetric solutions
Ψ0 = 1e-3*randn((n+1),(n+1))


# Get final solution
Ψ = CF.solve_timestepping(Ψ0, p, Δt, steps)
contourf(reverse(p.nodes), reverse(p.nodes), Ψ', aspect_ratio=1)



# Create gif of saved time steps 
# sol,time = CF.solve_timestepping_save(Ψ0, p, Δt, steps)

# plt = contourf(
#     xlim=(-1,1),
#     ylim=(-1,1),
#     title="Streamfunction Ψ on $n x $n",
#     titlefontsize=10,
#     xtickfontsize=8,
#     aspect_ratio=1,
#     legend=false,
# )
# plt2 = scatter(
#     xlim=(0,Δt*steps),
#     ylim=(-0.25,0.25),
#     title="Ψcenter value",
#     titlefontsize=10,
#     xtickfontsize=8,
#     ytickfontsize=8,
#     legend=false,
#     markercolor=:blue,
# )

# anim = @animate for i in 1:steps
#     plot(contourf!(plt, reverse(p.nodes), reverse(p.nodes), sol[i]'),
#     scatter!(plt2, [time[i]], [sol[i][nxc,nyc]], markercolor = :blue),
#         layout=grid(1,2, widths=(0.5,0.5)), size=(600,300)
#     )
# end
#
# gif(anim, "example_timestepping.gif", fps=15)
