println("--- Converge psis: ---")

# Create folder for converged solutions
folderconv_psis = "$(foldercont)/converged_psis"
mkdir(folderconv_psis)

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

# Stable asymmetric solutions
df_stab1 = filter(!(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn1.psi_c), df_pos)
df_stab2 = filter(!(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c), df_neg)

# Unstable asymmetric solutions
df_unstab1 = filter(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn1.psi_c, df)
df_unstab2 = filter(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c, df)


# Helper function to converge to exact Re values
# for asymmetric solutions
function converge_and_save_asymΨ(folder, folderconv_psis, name, Re, df, p)
    println("Converge asymmetric psi $name: Re = $Re")
    @unpack n,nodes = p.params

    # Find closest solution to use as initial guess
    i_cl = argmin(abs.(df.Re .- Re))
    row_cl = df[i_cl, :]

    p.params.Re = Re
    Ψ0 = readdlm("$folder/psis/psi_step$(@sprintf("%03d", row_cl.step))_Re$(@sprintf("%07.3f", row_cl.Re)).txt")

    Ψ, _, _ = CF.steadystate(Ψ0, p)
    writedlm("$folderconv_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)

    plt = contourf(reverse(nodes),
             reverse(nodes),
             Ψ',
             xlim = (-1, 1),
             ylim = (-1, 1),
             aspect_ratio = 1,
             axis = ([], false),
             legend = false,
             color = :davos,
             dpi = 800)
    fileplt = "$folderconv_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).png"
    savefig(plt, fileplt)
end

# Helper function to converge to symetric solutions
function converge_and_save_symΨ(folderconv_psis, name, Re, p)
    println("Converge symmetric psi $name: Re = $Re")
    @unpack n, nodes = p.params

    # Use zeros as intitial guess
    p.params.Re = Re
    Ψ0 = zeros(n + 1, n + 1)
    Ψ, _, _ = CF.steadystate(Ψ0, p)

    writedlm("$folderconv_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)

    # Save plot
    plt = contourf(reverse(nodes),
            reverse(nodes),
             Ψ',
             xlim = (-1, 1),
             ylim = (-1, 1),
             aspect_ratio = 1,
             axis = ([], false),
             legend = false,
             color = :davos,
             dpi = 800)
    fileplt = "$folderconv_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).png"
    savefig(plt, fileplt)
end

# Generate solutions
Re = 50
converge_and_save_symΨ(folderconv_psis, "sym", Re, p)

Re = 125
converge_and_save_asymΨ(foldercont, folderconv_psis, "asym1_stab", Re, df_stab1, p)
converge_and_save_asymΨ(foldercont, folderconv_psis, "asym2_stab", Re, df_stab2, p)
converge_and_save_symΨ(folderconv_psis, "sym", Re, p)

Re = 300
converge_and_save_asymΨ(foldercont, folderconv_psis, "asym1_stab", Re, df_stab1, p)
converge_and_save_asymΨ(foldercont, folderconv_psis, "asym2_stab", Re, df_stab2, p)
converge_and_save_asymΨ(foldercont, folderconv_psis, "asym1_unstab", Re, df_unstab1, p)
converge_and_save_asymΨ(foldercont, folderconv_psis, "asym2_unstab", Re, df_unstab2, p)
