using FourSidedCavityFlow
using CSV
using DataFrames
using DelimitedFiles
using Printf
using Plots
using LaTeXStrings
const CF = FourSidedCavityFlow
gr()

println("--- Branch 2 at saddle node: ---")

folderbranch2 = "$foldercont/branch2"
folderswitch = "$folderbranch2/switch_branch"
foldercont_branch2 = "$folderbranch2/continuation"
mkdir(folderbranch2)
mkdir(folderswitch)
mkdir(foldercont_branch2)

# Read results file of continuation as a DataFrame
df = CSV.read("$foldercont/results.csv", DataFrame)

# Saddle node is the starting point to switch branches
sn = df[argmax(df.Re), :]

# Compute an initial guess with continuation to start continuation 
# of branch 2
ΔRe = 1
steps = 6
Re_start = sn.Re

println("Arclength continuation to switch branch, ΔRe = $ΔRe, steps = $steps:")
Ψsn = readdlm("$foldercont/psis/psi_step$(@sprintf("%03d", sn.step))_Re$(@sprintf("%07.3f", sn.Re)).txt")
@time CF.continuation_arclength(folderswitch, Ψsn, p, Re_start, ΔRe, steps)

# Start continuation from the branch 2
df_switch = CSV.read("$(folderswitch)/results.csv", DataFrame)

ΔRe = -1
steps = 180
Re_start = 400
start = df_switch[argmin(abs.(df_switch.Re .- Re_start)), :]
folderbranch2 = "$foldercont/branch2"

println("Arclength continuation of branch 2, ΔRe = $ΔRe, steps = $steps:")
Ψi = readdlm("$folderswitch/psis/psi_step$(@sprintf("%03d", start.step))_Re$(@sprintf("%07.3f", start.Re)).txt")
@time CF.continuation_arclength(foldercont_branch2, Ψi, p, start.Re, ΔRe, steps)

# Plot
df_upper = filter(row -> row.Re > 200 && row.psi_c > 0, df)
df_branch2 = CSV.read("$foldercont_branch2/results.csv", DataFrame)

# Create 3D bifurcation curves
plt = plot(camera=(50,50),
           xtickfontsize = 7,
           ytickfontsize = 7,
           ztickfontsize = 7,
           xlabel = L"\mathrm{Re}",
           ylabel = L"u_{t}",
           zlabel = L"\Psi_{center}",
           xlabelfontsize = 12,
           ylabelfontsize = 12,
           zlabelfontsize = 12,
           thickness_scaling = 0.7,
           margin = 12Plots.mm,
           legend = false,
           dpi = 800)

plot!(df_branch2.Re, df_branch2.u_t, df_branch2.psi_c, linecolor = :green)
plot!(df_upper.Re, df_upper.u_t, df_upper.psi_c, linecolor = :black)

fileplt = "$foldercont/bifurcation_diag$(n)x$(n)_branch2.png"
savefig(plt, fileplt)
