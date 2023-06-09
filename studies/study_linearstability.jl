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
    @time CF.continuation_arclength_lsa(folderlsa, name, Ψi1, Ψi2, p, start1.Re, start2.Re, max_steps; Re_stop_mode = Re_stop_mode)
end

# Helper function to plot lambdas 
function plot_lambdas(df_lsa, folderlsa, name)
    c = palette(:Dark2_5).colors.colors

    plt = plot(xtickfontsize = 7,
                  ytickfontsize = 7,
                  thickness_scaling = 0.7,
                  xlabel = L"\mathrm{Re}",
                  ylabel = L"Real(\lambda)",
                  xlabelfontsize = 12,
                  ylabelfontsize = 12,
                  margin = 12Plots.mm,
                  legend = false,
                  palette = :Dark2_5, 
                  dpi = 800)

    plot!(plt, df_lsa.Re, df_lsa.lambda1re, marker = :circle, markersize = 3, linecolor = c[1], label = L"\lambda_1")
    plot!(plt, df_lsa.Re, df_lsa.lambda2re, marker = :cross, markersize = 4, linecolor = c[2], label = L"\lambda_2")
    plot!(plt, df_lsa.Re, df_lsa.lambda3re, marker = :x, markersize = 2, linecolor = c[3], label = L"\lambda_3")

    fileplt = "$folderlsa/lsa_$(name).png"
    savefig(plt, fileplt)
end

# Helper function to make a gif of real and imaginary part
# of lambdas
function gif_lambdas(folderlsa, name)
    # Open saved linearstability results
    df_lsa = CSV.read("$folderlsa/results_$(name).csv", DataFrame)

    c = palette(:Dark2_5).colors.colors
    ymax = maximum(maximum(eachcol(df_lsa[:, ["lambda1im", "lambda2im", "lambda3im"]])))

    if (ymax == 0)
        ymax = 0.01
    end

    plot(xlim = (-0.1, 0.1),
         ylim = (-ymax-0.2*ymax, ymax+0.2*ymax),
                  ytickfontsize = 7,
                  thickness_scaling = 0.7,
                  xlabel = L"Real(\lambda)",
                  ylabel = L"Imag(\lambda)",
                  xlabelfontsize = 12,
                  ylabelfontsize = 12,
                  margin = 12Plots.mm,
                  legend = false,
                  dpi = 800)

    anim = @animate for i in 1:size(df_lsa,1)
        scatter!([df_lsa.lambda1re[i]], [df_lsa.lambda1im[i]], color = c[1])
        scatter!([df_lsa.lambda2re[i]], [df_lsa.lambda2im[i]], color = c[2])
        scatter!([df_lsa.lambda3re[i]], [df_lsa.lambda3im[i]], color = c[3])
    end

    filegif = "$folderlsa/lsa_$(name).gif"
    gif(anim, filegif, fps=6)
end

# Helper function to start LSA
function get_Re_start(df, Re_start; incr = true)
    i_start1 = argmin(abs.(df.Re .- Re_start))
    start1 = df[i_start1, :]
    start2a = filter(row -> row.step == (start1.step + 1), df)[1, :]
    start2b = filter(row -> row.step == (start1.step - 1), df)[1, :]

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
max_steps = 50

### Saddle node ###
name = "sn"
# Re_start = 346 
# start1, start2 = get_Re_start(df, Re_start)
# lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps, 1, p)

# Plot approaching from above and below
df_lsa = CSV.read("$folderlsa/results_$(name).csv", DataFrame)
i_sn = argmax(df_lsa.Re)
plot_lambdas(df_lsa[1:i_sn, :], folderlsa, "$(name)_1")
plot_lambdas(df_lsa[i_sn:end, :], folderlsa, "$(name)_2")
gif_lambdas(folderlsa, name)

### Pitchfork 1 ###
name = "pf1"
# df_1 = filter(row -> row.Re < 100, df)
# Re_start = 67
# start1, start2 = get_Re_start(df_1, Re_start; incr = false)
# lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps, 2, p)

df_lsa = CSV.read("$folderlsa/results_$(name).csv", DataFrame)
i_pf1 = argmin(df_lsa.Re)
df_lsa1 = df_lsa[1:i_pf1, :]
plot_lambdas(df_lsa1, folderlsa, name)
gif_lambdas(folderlsa, name)

### Pitchfork 2 ###
name = "pf2"
# df_2 = filter(row -> row.Re > 100, df)
# Re_start = 173.5
# start1, start2 = get_Re_start(df_2, Re_start; incr = false)
# lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps, 2, p)

df_lsa = CSV.read("$folderlsa/results_$(name).csv", DataFrame)
i_pf2 = argmin(df_lsa.Re)
plot_lambdas(df_lsa[1:i_pf2, :], folderlsa, name)
gif_lambdas(folderlsa, name)

### Branch 2 ###
name = "sn_branch2"
# df_branch2 = CSV.read("$foldercont_branch2/results.csv", DataFrame)
# Re_start = 354
# start1, start2 = get_Re_start(df_branch2, Re_start; incr = false)
# lsa_around_bif_point(foldercont_branch2, folderlsa, name, start1, start2, max_steps, 2, p)

df_lsa = CSV.read("$folderlsa/results_$(name).csv", DataFrame)
i_sn = argmin(df_lsa.Re)
plot_lambdas(df_lsa[1:i_sn, :], folderlsa, name)
gif_lambdas(folderlsa, name)
