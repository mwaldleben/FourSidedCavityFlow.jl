# Helper function to order the results of the continuation algorithm
function order_cont_results(foldercont)
    # Read results file of continuation as a DataFrame
    df = CSV.read("$foldercont/results.csv", DataFrame)

    # The DataFrame has to be filtered to distinguish
    # between stable and unstable solutions
    # u : upper solution, l : lower solution
    df_u = filter(row -> row.psi_c > 0, df)
    df_l = filter(row -> row.psi_c < 0, df)

    # The approximated saddle node bifurcations have maximum Re values
    sn_u = df_u[argmax(df_u.Re), :]
    sn_l = df_l[argmax(df_l.Re), :]

    # Approximated pitchfork 1
    df1 = filter(row -> row.Re < 100, df)
    i_pf1 = argmin(abs.(df1.psi_c))
    pf1 = df1[i_pf1, :]

    # Approximated pitchfork 2
    df2 = filter(row -> row.Re > 100, df)
    i_pf2 = argmin(abs.(df2.psi_c))
    pf2 = df2[i_pf2, :]

    # Roughly where hopf bifurcation is 
    Re_hopf = 348

    # Unstable asymmetric solutions
    df_unstab_u = filter(row -> row.Re > 150 && pf2.psi_c <= row.psi_c <= sn_u.psi_c, df)
    df_unstab_l = filter(row -> row.Re > 150 && sn_l.psi_c <= row.psi_c <= pf2.psi_c, df)

    df_unstab_hopf_u = filter(row -> row.Re >= Re_hopf && row.psi_c > sn_u.psi_c, df)
    df_unstab_hopf_l = filter(row -> row.Re >= Re_hopf && row.psi_c < sn_l.psi_c, df)

    # Stable asymmetric solutions
    df_stab_ua = filter(row -> 150 <= row.Re <= Re_hopf && row.psi_c > sn_u.psi_c, df_u)
    df_stab_ub = filter(row -> row.Re <= 160, df_u)
    df_stab_u = unique(vcat(df_stab_ua, df_stab_ub))
    sort!(df_stab_u, "Re")

    df_stab_la = filter(row -> 150 <= row.Re <= Re_hopf && row.psi_c < sn_l.psi_c, df_l)
    df_stab_lb = filter(row -> row.Re <= 160, df_l)
    df_stab_l = unique(vcat(df_stab_la, df_stab_lb))
    sort!(df_stab_l, "Re")

    return df,
    df_u,
    df_l,
    df_unstab_u,
    df_unstab_l,
    df_unstab_hopf_u,
    df_unstab_hopf_l,
    df_stab_u,
    df_stab_l,
    sn_u,
    sn_l,
    pf1,
    pf2
end

# Helper function to converge to exact Re values
# for asymmetric solutions
function converge_and_save_asymΨ(foldercont, folderpsis, name, Re, df, p)
    println("Converge asymmetric psi $name: Re = $Re")
    @unpack n, nodes = p.params

    # Find closest solution to use as initial guess
    i_cl = argmin(abs.(df.Re .- Re))
    row_cl = df[i_cl, :]

    p.params.Re = Re
    Ψ0 = readdlm("$foldercont/psis/psi_step$(@sprintf("%03d", row_cl.step))_Re$(@sprintf("%07.3f", row_cl.Re)).txt")

    Ψ, _, _ = CF.steadystate(Ψ0, p)
    writedlm("$folderpsis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)

    plt = contourf(reverse(nodes),
        reverse(nodes),
        Ψ';
        xlim = (-1, 1),
        ylim = (-1, 1),
        aspect_ratio = 1,
        axis = ([], false),
        legend = false,
        color = :davos,
        dpi = 800)
    fileplt = "$folderpsis/psi_Re$(@sprintf("%07.3f", Re))_$(name).png"
    savefig(plt, fileplt)

    return nothing
end

# Helper function to converge to symetric solutions
function converge_and_save_symΨ(folderpsis, name, Re, p)
    println("Converge symmetric psi $name: Re = $Re")
    @unpack n, nodes = p.params

    # Use zeros as intitial guess
    p.params.Re = Re
    Ψ0 = zeros(n + 1, n + 1)
    Ψ, _, _ = CF.steadystate(Ψ0, p)

    writedlm("$folderpsis/psi_Re$(@sprintf("%07.3f", Re))_$(name).txt", Ψ)

    # Save plot
    plt = contourf(reverse(nodes),
        reverse(nodes),
        Ψ';
        xlim = (-1, 1),
        ylim = (-1, 1),
        aspect_ratio = 1,
        axis = ([], false),
        legend = false,
        color = :davos,
        dpi = 300)
    fileplt = "$folderpsis/psi_Re$(@sprintf("%07.3f", Re))_$(name).png"
    savefig(plt, fileplt)
    
    return nothing
end

# Helper function to do stability analysis around bifurcation point
function lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps,
    Re_stop_mode, p)
    Ψi1 = readdlm("$foldercont/psis/psi_step$(@sprintf("%03d", start1.step))_Re$(@sprintf("%07.3f", start1.Re)).txt")
    Ψi2 = readdlm("$foldercont/psis/psi_step$(@sprintf("%03d", start2.step))_Re$(@sprintf("%07.3f", start2.Re)).txt")

    println("Arclength continuation with linear stability analysis $(name), max_steps = $max_steps:")
    @time CF.continuation_arclength_lsa(folderlsa,
        name,
        Ψi1,
        Ψi2,
        p,
        start1.Re,
        start2.Re,
        max_steps;
        Re_stop_mode = Re_stop_mode)

    return nothing
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
