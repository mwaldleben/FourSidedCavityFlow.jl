import FourSidedCavityFlow as CF
using DelimitedFiles
using Random
using Plots
gr()

# Read a steady state Ψ generated with "study_steadystate.jl"
n = 32
Re = 320
p = CF.setup_params(n, Re)
name = "asym_upper1"

title = "--- Linear stability analysis: $(n)x$(n) ---"
println(title)

filename = "results/$(n)x$(n)_$(name)_Re$(Re).txt"
filepath = joinpath(@__DIR__, filename)
Ψstart = readdlm(filepath)
println("Reading Ψstart from $(filename)")

Re_start = Re
ΔRe = +1
steps = 350

nxc = Int(floor(n / 2) + 1)
nyc = Int(floor(n / 2) + 1)

println("Continuation, $(steps) steps:")
sol, lambdas, Re_series = @time CF.continuation_arclength(Ψstart, p, Re_start, ΔRe, steps;
                                                          linearstability = true)

Ψcenter = Vector{Float64}(undef, size(sol)[1])
for i in 1:size(sol)[1]
    Ψcenter[i] = sol[i][nxc, nyc]
end

plt = plot(Re_series,
           Ψcenter,
           xlim = (0, 400),
           ylim = (-0.25, 0.25),
           titlefontsize = 10,
           xtickfontsize = 9,
           ytickfontsize = 9,
           xlabel = "Re",
           ylabel = "Ψ center value",
           xlabelfontsize = 10,
           ylabelfontsize = 10,
           margin = 12Plots.mm,
           legend = false,
           palette = :davos,
           dpi = 600)

pltname = "results/$(n)x$(n)_linearstability_p2_1.png"
pltpath = joinpath(@__DIR__, pltname)
savefig(plt, pltpath)
println("Saved plot 1 to $(pltname)")

lambda1 = similar(Re_series)
lambda2 = similar(Re_series)
lambda3 = similar(Re_series)
for i in 1:steps
    lambda1[i] = lambdas[i][1]
    lambda2[i] = lambdas[i][2]
    lambda3[i] = lambdas[i][3]
end

plt = scatter(Re_series,
              lambda1,
              xlim = (200, 400),
              titlefontsize = 10,
              xtickfontsize = 9,
              ytickfontsize = 9,
              xlabel = "Re",
              ylabel = "Real(λ)",
              xlabelfontsize = 10,
              ylabelfontsize = 10,
              margin = 20Plots.mm,
              markershape = :circle,
              markersize = 2,
              palette = :davos,
              label = "λ1",
              dpi = 600)

scatter!(Re_series, lambda2, marker = :cross, markersize = 3, label = "λ2")
scatter!(Re_series, lambda3, marker = :star5, markersize = 3, label = "λ2 ")

pltname = "results/$(n)x$(n)_linearstability_p2_2.png"
pltpath = joinpath(@__DIR__, pltname)
savefig(plt, pltpath)
println("Saved plot 2, to $(pltname)")
println(repeat('-', length(title)))
