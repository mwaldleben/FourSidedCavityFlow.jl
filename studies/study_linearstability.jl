println("--- Linear stability analysis: ---")

foldercont_branch2 = "$foldercont/branch2/continuation"
folderlsa = "$(foldercont)/lsa"
# mkdir(folderlsa)

# Read results file of continuation as a DataFrame
df = CSV.read("$foldercont/results.csv", DataFrame)

df_pos = filter(row -> row.psi_c > 0, df)
df_neg = filter(row -> row.psi_c < 0, df)

# Helper function to do stability analysis around bifurcation point
function lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps, Re_stop_mode, p)
    Ψi1 = readdlm("$foldercont/psis/psi_step$(@sprintf("%03d", start1.step))_Re$(@sprintf("%07.3f", start1.Re)).txt")
    Ψi2 = readdlm("$foldercont/psis/psi_step$(@sprintf("%03d", start2.step))_Re$(@sprintf("%07.3f", start2.Re)).txt")

    println("Arclength continuation with linear stability analysis $(name), max_steps = $max_steps:")
    CF.continuation_arclength_lsa(folderlsa, name, Ψi1, Ψi2, p, start1.Re, start2.Re, max_steps; Re_stop_mode = Re_stop_mode)
end

# Helper function to plot lambdas 
function plot_lambdas(foldercont, folderlsa, name)
    # Open saved linearstability results
    df_lsa = CSV.read("$folderlsa/results_$(name).csv", DataFrame)

    plt = plot(xtickfontsize = 7,
                  ytickfontsize = 7,
                  thickness_scaling = 0.7,
                  xlabel = L"\mathrm{Re}",
                  ylabel = L"Real(\lambda)",
                  xlabelfontsize = 12,
                  ylabelfontsize = 12,
                  margin = 12Plots.mm,
                  legend = false,
                  markershape = :circle,
                  markersize = 3,
                  palette = :Dark2_5, 
                  dpi = 800)

    plot!(df_lsa.Re, df_lsa.lambda1, marker = :circle, markersize = 3, label = L"\lambda_1")
    plot!(df_lsa.Re, df_lsa.lambda2, marker = :cross, markersize = 4, label = L"\lambda_2")
    plot!(df_lsa.Re, df_lsa.lambda3, marker = :x, markersize = 2, label = L"\lambda_3")

    fileplt = "$folderlsa/lsa_$(name).png"
    savefig(plt, fileplt)
end

# Helper function to start LSA
function get_Re_start(df, bif, Re_start; side = true, upper = true, incr = true)
    if side == true 
        if upper == true 
            df_filt = filter(row -> row.psi_c > bif.psi_c, df)
        else
            df_filt = filter(row -> row.psi_c < bif.psi_c, df)
        end
    else
        df_filt = df
    end

    i_start1 = argmin(abs.(df_filt.Re .- Re_start))
    start1 = df_filt[i_start1, :]
    start2a = filter(row -> row.step == (start1.step + 1), df_filt)[1, :]
    start2b = filter(row -> row.step == (start1.step - 1), df_filt)[1, :]

    if start2a.Re > start1.Re
        if incr == true
            start2 = start2a
        else
            start2 = start2b
        end
    else
        if incr == true
            start2 = start2b
        else
            start2 = start2a
        end
    end

    return start1, start2
end

# Linear stability analysis
# Set max steps for lsa 
max_steps = 10

### Saddle node ###
sn = df_pos[argmax(df_pos.Re), :]

# Approaching from above (upper)
# Re_start_u = 346 
# start1, start2 = get_Re_start(df, sn, Re_start_u)
# lsa_around_bif_point(foldercont, folderlsa, "sn_u", start1, start2, max_steps, 1, p)
# plot_lambdas(foldercont, folderlsa, "sn_u")

# Approaching from below (lower)
# Re_start_u = 351 
# start1, start2 = get_Re_start(df, sn, Re_start_u; upper = false)
# lsa_around_bif_point(folderlsa, foldercont, "sn_l", start1, start2, max_steps, 1, p)
# plot_lambdas(foldercont, folderlsa, "sn_l")

### Pitchfork 1 ###
# df_1 = filter(row -> row.Re < 100, df)
# pf1 = df_1[argmin(abs.(df_1.Re)), :]
# Re_start_u = 66.4
# start1, start2 = get_Re_start(df_1, pf1, Re_start_u; incr = false)
# lsa_around_bif_point(foldercont, foldercont, "pf1", start1, start2, max_steps, 2, p)
# plot_lambdas(foldercont, folderlsa, "pf1")

### Pitchfork 2 ###
# df_2 = filter(row -> row.Re > 100, df)
# pf2 = df_2[argmin(abs.(df_2.Re)), :]
# Re_start_u = 172.8
# start1, start2 = get_Re_start(df_2, pf2, Re_start_u; incr = false)
# lsa_around_bif_point(foldercont, foldercont, "pf2", start1, start2, max_steps, 2, p)
# plot_lambdas(foldercont, folderlsa, "pf2")

### Branch 2 ###
df_branch2 = CSV.read("$foldercont_branch2/results.csv", DataFrame)
Re_start = 352
sn_branch2 = df_branch2[argmin(df_branch2.Re), :]
start1, start2 = get_Re_start(df_branch2, sn_branch2, Re_start; side = false, incr = false)
lsa_around_bif_point(foldercont_branch2, folderlsa, "sn_branch2", start1, start2, max_steps, 0, p)
plot_lambdas(foldercont, folderlsa, "sn_branch2")
