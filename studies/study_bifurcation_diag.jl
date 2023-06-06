using FourSidedCavityFlow
using Plots
using CSV
using DataFrames
using Random
using LaTeXStrings
using UnPack
using Printf
using DelimitedFiles
const CF = FourSidedCavityFlow
gr()

# Read results file of continuation as a DataFrame
foldername = "continuation$(n)x$(n)"

df = CSV.read("$foldername/results.csv", DataFrame)

# The DataFrame has to be filtered to create
# the bifurcation diagram
df_pos = filter(row -> row.psi_c > 0, df)
df_neg = filter(row -> row.psi_c < 0, df)

# The saddle node bifurcations have maximum Re values
sn1 = df_pos[argmax(df_pos.Re), :]
sn2 = df_neg[argmax(df_neg.Re), :]

# Pitchfork 1
df_1 = filter(row -> row.Re < 100, df)
row_pf1 = argmin(abs.(df_1[!, "psi_c"]))
pf1 = df_1[row_pf1, :]

# Pitchfork 2
df_2 = filter(row -> row.Re > 100, df)
row_pf2 = argmin(abs.(df_2[!, "psi_c"]))
pf2 = df_2[row_pf2, :]

# Create bifurcation curve
plt = plot(xlim = (0, 400),
           ylim = (-0.35, 0.35),
           xtickfontsize = 7,
           ytickfontsize = 7,
           xlabel = L"\mathrm{Re}",
           ylabel = L"\Psi_{center}",
           xlabelfontsize = 12,
           ylabelfontsize = 12,
           thickness_scaling = 0.6,
           margin = 12Plots.mm,
           legend = false,
           grid=false,
           palette = :davos,
           dpi = 600)

# Plot stable assymetric solutions
df_stab1 = filter(!(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn1.psi_c), df_pos)
sort!(df_stab1, ["Re"])
plot!(df_stab1[!, "Re"], df_stab1[!, "psi_c"])

df_stab2 = filter(!(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c), df_neg)
sort!(df_stab2, ["Re"])
plot!(df_stab2[!, "Re"], df_stab2[!, "psi_c"])

# Plot unstable assymetric solutions
df_unstab1 = filter(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn1.psi_c, df)
sort!(df_unstab1, ["Re"])
plot!(df_unstab1[!, "Re"], df_unstab1[!, "psi_c"], linestyle = :dot)

df_unstab2 = filter(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c, df)
sort!(df_unstab2, ["Re"])
plot!(df_unstab2[!, "Re"], df_unstab2[!, "psi_c"], linestyle = :dot)

# Plot stable symmetric solutions
Re_stab1 = range(0, pf1.Re, length = 100)
ψstab1 = zeros(size(Re_stab1))
plot!(Re_stab1, ψstab1, linewidth = 1, palette = :davos)

Re_stab2 = range(pf2.Re, 400, length = 100)
ψstab2 = zeros(size(Re_stab2))
plot!(Re_stab2, ψstab2, palette = :davos)

# Plot unstable symmetric solution
Re_unstab = range(pf1.Re, pf2.Re, length = 100)
ψunstab = zeros(size(Re_unstab))
plot!(Re_unstab, ψunstab, linestyle = :dot, palette = :davos)


# Add solutions as inset plots (optional)
function inset_plot(foldername, name, Re, p, pos_relx, pos_rely, subplot_nb)
    @unpack nodes, ic = p.params

    Ψ = readdlm("$(foldername)/converged_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt")

    contour!(reverse(nodes),
             reverse(nodes),
             Ψ',
             # color=[:black, :black],
             color=:davos,
             title = L"\mathrm{Re} \; %$(Re)",
             titlefontsize = 9,
             inset = (1, bbox(pos_relx, pos_rely, 0.13, 0.13, :bottom, :left)),
             xlim = (-1, 1),
             ylim = (-1, 1),
             aspect_ratio = 1,
             axis = ([], false),
             legend = false,
             subplot = subplot_nb,
             bg_inside = "#f2f2f2",
            )
    scatter!([Re], [Ψ[ic, ic]], marker = :circle, markersize = 2, palette = :davos)
end

n = 32
p = CF.setup_struct(n, 0)


Re = 50
inset_plot(foldername, "sym", Re, p, 0.05, 0.52, 2)

Re = 125
inset_plot(foldername, "sym", Re, p, 0.255, 0.52, 3)
inset_plot(foldername, "asym1_stab", Re, p, 0.255, 0.83, 4)
inset_plot(foldername, "asym2_stab", Re, p, 0.255, 0.025, 5)

Re = 300
inset_plot(foldername, "asym1_stab", Re, p, 0.685, 0.8, 6)
inset_plot(foldername, "asym2_stab", Re, p, 0.685, 0.035, 7)
inset_plot(foldername, "asym1_unstab", Re, p, 0.685, 0.533, 8)
inset_plot(foldername, "asym2_unstab", Re, p, 0.685, 0.335, 9)

# Save plot
fileplt = "$foldername/bifurcation_diag$(n)x$(n).png"
savefig(plt, fileplt)
println("Saved plot to $(fileplt)")
gui()
