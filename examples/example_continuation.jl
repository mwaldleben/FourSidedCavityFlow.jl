import CavityFlow as CF
using Random
using Plots
gr()

Random.seed!(1234)

n = 16
Re = 100
p = CF.setup_params(n, Re)

Δt = 1
time_steps = 140

# Convergence to one of the two possible assymetric solutions
Ψ0 = 1e-3 * randn((n + 1), (n + 1))

# Find an unstable solution
Ψ0 = CF.solve_timestepping(Ψ0, p, Δt, time_steps)

Re_start = Re
ΔRe = -1
steps = 800

nxc = Int(floor(n / 2) + 1)
nyc = Int(floor(n / 2) + 1)

sol, Re_series = @time CF.solve_continuation(Ψ0, p, Re_start, ΔRe, steps)

plt = contourf(xlim = (-1, 1),
               ylim = (-1, 1),
               title = "Streamfunction Ψ on $n x $n",
               titlefontsize = 10,
               xtickfontsize = 8,
               aspect_ratio = 1,
               legend = false)
plt2 = scatter(xlim = (0, 300),
               ylim = (-0.25, 0.25),
               title = "Ψcenter value",
               titlefontsize = 10,
               xtickfontsize = 8,
               ytickfontsize = 8,
               legend = false,
               markercolor = :blue)

@time anim = @animate for i in 1:steps
    println(i)
    plot(contourf(plt, reverse(p.nodes), reverse(p.nodes), sol[i]'),
         scatter!(plt2, [Re_series[i]], [sol[i][nxc, nyc]], markercolor = :blue),
         layout = grid(1, 2, widths = (0.5, 0.5)), size = (600, 300))
end
gif(anim, "example_continuation.gif", fps = 45)
