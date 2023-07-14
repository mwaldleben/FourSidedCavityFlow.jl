println("--- Study branch 2 at saddle node: ---")

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

### Linearstability Branch 2 ### 
name = "branch2"
max_steps = 100
df_branch2 = CSV.read("$foldercont_branch2/results.csv", DataFrame)
Re_start = 354
start1, start2 = get_Re_start(df_branch2, Re_start; incr = false)
lsa_around_bif_point(foldercont_branch2, folderlsa, name, start1, start2, max_steps, 2, p)
