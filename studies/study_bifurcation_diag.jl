println("--- Generate bifurcation diagram: ---")

# Read results file of continuation as a DataFrame
df = CSV.read("$foldercont/results.csv", DataFrame)

# The DataFrame has to be filtered to distinguish
# between stable and unstable solutions
df_pos = filter(row -> row.psi_c > 0, df)
df_neg = filter(row -> row.psi_c < 0, df)

# The approximated saddle node bifurcations have maximum Re values
sn1 = df_pos[argmax(df_pos.Re), :]
sn2 = df_neg[argmax(df_neg.Re), :]

# Approximated pitchfork 1
df1 = filter(row -> row.Re < 100, df)
i_pf1 = argmin(abs.(df1.psi_c))
pf1 = df1[i_pf1, :]

# Approximated pitchfork 2
df2 = filter(row -> row.Re > 100, df)
i_pf2 = argmin(abs.(df2.psi_c))
pf2 = df2[i_pf2, :]

# Create bifurcation curve
plt = plot(;
    xlim = (0, 400),
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
    grid = false,
    dpi = 800)

# Plot unstable asymmetric solutions
Re_hopf = 348 # TODO: has to be calculated
psi_c_hopf = 0.1690
df_unstab1 = filter(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn1.psi_c, df)
plot!(df_unstab1.Re, df_unstab1.psi_c; linestyle = :dot, linecolor = :black)

df_unstab2 = filter(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c, df)
plot!(df_unstab2.Re, df_unstab2.psi_c; linestyle = :dot, linecolor = :black)

df_unstab_hopf1 = filter(row -> row.Re >= Re_hopf && row.psi_c > sn1.psi_c, df)
plot!(df_unstab_hopf1.Re, df_unstab_hopf1.psi_c; linestyle = :dot, linecolor = :black)

df_unstab_hopf2 = filter(row -> row.Re >= Re_hopf && row.psi_c < sn2.psi_c, df)
plot!(df_unstab_hopf2.Re, df_unstab_hopf2.psi_c; linestyle = :dot, linecolor = :black)

df_unstab_hopf2 = filter(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c, df)
plot!(df_unstab2.Re, df_unstab2.psi_c; linestyle = :dot, linecolor = :black)

# Plot stable asymmetric solutions
# sort!(df_stab1, "Re")
df_stab1a = filter(row -> 150 <= row.Re <= Re_hopf && row.psi_c > sn1.psi_c, df_pos)
df_stab1b = filter(row -> row.Re <= 160, df_pos)
sort!(df_stab1a, "Re")
sort!(df_stab1b, "Re")
plot!(df_stab1a.Re, df_stab1a.psi_c; linecolor = :black)
plot!(df_stab1b.Re, df_stab1b.psi_c; linecolor = :black)

df_stab2a = filter(row -> 150 <= row.Re <= Re_hopf && row.psi_c < sn2.psi_c, df_neg)
df_stab2b = filter(row -> row.Re <= 160, df_neg)
sort!(df_stab2a, "Re")
sort!(df_stab2b, "Re")
plot!(df_stab2a.Re, df_stab2a.psi_c; linecolor = :black)
plot!(df_stab2b.Re, df_stab2b.psi_c; linecolor = :black)

# Plot stable symmetric solutions
Re_stab1 = range(0, pf1.Re; length = 100)
ψstab1 = zeros(size(Re_stab1))
plot!(Re_stab1, ψstab1; linecolor = :black)

Re_stab2 = range(pf2.Re, 400; length = 100)
ψstab2 = zeros(size(Re_stab2))
plot!(Re_stab2, ψstab2; linecolor = :black)

# Plot unstable symmetric solution
Re_unstab = range(pf1.Re, pf2.Re; length = 100)
ψunstab = zeros(size(Re_unstab))
plot!(Re_unstab, ψunstab; linestyle = :dot, linecolor = :black)

# Helper function to plot solutions as small inset frames (optional)
function inset_plot(folder, name, Re, p, pos_relx, pos_rely, subplot_nb)
    @unpack nodes, ic = p.params

    Ψ = readdlm("$folder/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt")

    contour!(reverse(nodes),
        reverse(nodes),
        Ψ';
        title = L"\mathrm{Re} \; %$(Re)",
        titlefontsize = 9,
        inset = (1, bbox(pos_relx, pos_rely, 0.13, 0.13, :bottom, :left)),
        xlim = (-1, 1),
        ylim = (-1, 1),
        aspect_ratio = 1,
        axis = ([], false),
        legend = false,
        subplot = subplot_nb,
        color = :davos,
        bg_inside = "#f2f2f2",
        dpi = 800)
    return scatter!([Re], [Ψ[ic, ic]]; marker = :circle, markersize = 2, color = :black)
end

Re = 50
inset_plot(folderconv_psis, "sym", Re, p, 0.05, 0.52, 2)

Re = 125
inset_plot(folderconv_psis, "sym", Re, p, 0.255, 0.52, 3)
inset_plot(folderconv_psis, "asym1_stab", Re, p, 0.255, 0.83, 4)
inset_plot(folderconv_psis, "asym2_stab", Re, p, 0.255, 0.025, 5)

Re = 300
inset_plot(folderconv_psis, "asym1_stab", Re, p, 0.685, 0.8, 6)
inset_plot(folderconv_psis, "asym2_stab", Re, p, 0.685, 0.035, 7)
inset_plot(folderconv_psis, "asym1_unstab", Re, p, 0.685, 0.533, 8)
inset_plot(folderconv_psis, "asym2_unstab", Re, p, 0.685, 0.335, 9)

# Save plot
fileplt = "$foldercont/bifurcation_diag$(n)x$(n).png"
savefig(plt, fileplt)
