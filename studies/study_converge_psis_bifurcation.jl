using FourSidedCavityFlow
using Plots
using CSV
using DataFrames
using LaTeXStrings
using Printf
using DelimitedFiles
const CF = FourSidedCavityFlow
gr()

# Read results file of continuation as a DataFrame
foldername = "continuation$(n)x$(n)"

df = CSV.read("$foldername/results.csv", DataFrame)

# The DataFrame has to be filtered to find
# the approximative bifurcation points
df_pos = filter(row -> row.psi_c > 0, df)
df_neg = filter(row -> row.psi_c < 0, df)

# The approximated saddle node bifurcations have maximum Re values
sn1 = df_pos[argmax(df_pos.Re), :]
sn2 = df_neg[argmax(df_neg.Re), :]

# Approximated Pitchfork 1
df_1 = filter(row -> row.Re < 100, df)
row_pf1 = argmin(abs.(df_1[!, "psi_c"]))
pf1 = df_1[row_pf1, :]

# Approximated Pitchfork 2
df_2 = filter(row -> row.Re > 100, df)
row_pf2 = argmin(abs.(df_2[!, "psi_c"]))
pf2 = df_2[row_pf2, :]

# Converge bifurcation points
# Helper function to converge Ψ with 1D newton
function converge_and_saveΨ_bif(foldername, name, df_row, p)
    println("Converge bifurcation point $name : Re_start = $(df_row.Re)")
    Ψ0 = readdlm("$(foldername)/psis/psi_step$(@sprintf("%03d", df_row.step))_Re$(@sprintf("%07.3f", df_row.Re)).txt")
    u0 = Ψ0[3:(n - 1), 3:(n - 1)][:]

    Re0 = df_row.Re
    Re, u, iter, tol = CF.newton1D_for_linearstability(Re0, u0, p; tolmax = 1e-12,
                                                       maxiter = 20)
    p.params.Re = Re
    Ψ = CF.constructBC(u, p)

    println("  $name : Re = $Re")
    writedlm("$(foldername)/converged_psis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)
end

# Run 1D newton
converge_and_saveΨ_bif(foldername, "pf1", pf1, p)
converge_and_saveΨ_bif(foldername, "pf2", pf2, p)
