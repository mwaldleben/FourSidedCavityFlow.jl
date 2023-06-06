using FourSidedCavityFlow
using CSV
using DataFrames
using Printf
using UnPack
using DelimitedFiles
const CF = FourSidedCavityFlow

# Read results file of continuation as a DataFrame
# and create folder for converged solutions
foldername = "continuation$(n)x$(n)"
folderconv_psis = "$(foldername)/converged_psis"
mkdir(folderconv_psis)

df = CSV.read("$foldername/results.csv", DataFrame)

# The DataFrame has to be filtered to create 
# to distinguish between stable and unstable solutions
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

# stable assymetric solutions
df_stab1 = filter(!(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn1.psi_c), df_pos)
df_stab2 = filter(!(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c), df_neg)

# unstable assymetric solutions
df_unstab1 = filter(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn1.psi_c, df)
df_unstab2 = filter(row -> row.Re > 150 && sn2.psi_c <= row.psi_c <= pf2.psi_c, df)


# Helper function to converge to exact Re values
# for asymmetric solutions
function converge_and_save_asymΨ(foldername, name, Re, df, p)
    @unpack n = p.params

    # Find closest solution to use as initial guess
    closest_row = argmin(abs.(df.Re .- Re))
    df_row = df[closest_row, :]

    p.params.Re = Re
    Ψ0 = readdlm("$(foldername)/psis/psi_step$(@sprintf("%03d", df_row.step))_Re$(@sprintf("%07.3f", df_row.Re)).txt")

    Ψ, iter, tol = CF.steadystate(Ψ0, p)

    writedlm("$(foldername)/converged_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)
end

# Helper function to converge to symetric solutions
function converge_and_save_symΨ(foldername, name, Re, p)
    @unpack n = p.params

    p.params.Re = Re
    Ψ0 = zeros(n + 1, n + 1)
    Ψ, iter, tol = CF.steadystate(Ψ0, p)

    writedlm("$(foldername)/converged_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)
end

# Generate converged solutions
Re = 50
converge_and_save_symΨ(foldername, "sym", Re, p)

Re = 125
converge_and_save_asymΨ(foldername, "asym1_stab", Re, df_stab1, p)
converge_and_save_asymΨ(foldername, "asym2_stab", Re, df_stab2, p)
converge_and_save_symΨ(foldername, "sym", Re, p)

Re = 300
converge_and_save_asymΨ(foldername, "asym1_stab", Re, df_stab1, p)
converge_and_save_asymΨ(foldername, "asym2_stab", Re, df_stab2, p)
converge_and_save_asymΨ(foldername, "asym1_unstab", Re, df_unstab1, p)
converge_and_save_asymΨ(foldername, "asym2_unstab", Re, df_unstab2, p)
