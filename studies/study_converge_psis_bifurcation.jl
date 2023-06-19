println("--- Converge psis of bifurcations: ---")

# Folder for converged solutions
folderconv_psis = "$(foldercont)/converged_psis"

# Read results file of continuation as a DataFrame
df = CSV.read("$foldercont/results.csv", DataFrame)

# The DataFrame has to be filtered to distinguish
# between stable and unstable solutions
df_pos = filter(row -> row.psi_c > 0, df)
df_neg = filter(row -> row.psi_c < 0, df)

# The approximated saddle node bifurcations have maximum Re values
sn = df[argmax(df_pos.Re), :]

# Approximated pitchfork 1
df1 = filter(row -> row.Re < 100, df)
i_pf1 = argmin(abs.(df1.psi_c))
pf1 = df1[i_pf1, :]

# Approximated pitchfork 2
df2 = filter(row -> row.Re > 100, df)
i_pf2 = argmin(abs.(df2.psi_c))
pf2 = df2[i_pf2, :]

# Helper function to converge to Ψ having eigenvalue 0
# with 1D newton
function converge_and_saveΨ_bif(folder, folderconv_psis, name, row_bif, p)
    println("Converge bifurcation point $name: Re_start = $(row_bif.Re)")
    @unpack nodes = p.params

    Ψ0 = readdlm("$folder/psis/psi_step$(@sprintf("%03d", row_bif.step))_Re$(@sprintf("%07.3f", row_bif.Re)).txt")
    u0 = Ψ0[3:(n - 1), 3:(n - 1)][:]

    Re0 = row_bif.Re
    Re, u, _, _ = CF.newton1D_for_linearstability(Re0, u0, p; tolmax = 1e-12, maxiter = 10,
        verbose = true)
    p.params.Re = Re
    Ψ = CF.constructBC(u, p)

    println("  $name : Re = $Re")
    return writedlm("$folderconv_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)
end

# Converge bifurcation points
converge_and_saveΨ_bif(foldercont, folderconv_psis, "pf1", pf1, p)
converge_and_saveΨ_bif(foldercont, folderconv_psis, "pf2", pf2, p)
converge_and_saveΨ_bif(foldercont, folderconv_psis, "sn", sn, p)
