function continuation_arclength(Ψi, p::CavityStruct, Re_start, ΔRe, steps; save_steps = 1)
    @unpack n, scl = p.params

    # Create directories
    foldername = "continuation$(n)x$(n)"
    folderpsis = "$(foldername)/psis"
    mkdir(foldername)
    mkdir(folderpsis)

    # Write header and open results file
    fileresults = "$(foldername)/results.csv"
    header = ["step", "Re", "psi_c", "psi_ul", "psi_ur", "psi_ll", "psi_lr", "newton_iterations", "newton_time"]
    writedlm(fileresults, reshape(header,1,length(header)), ',')
    io = open(fileresults, "a")

    p.params.Re = Re_start
    @inbounds u0 = reshape(Ψi[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    u1, _, _ = newton(f!, u0, p)
    Ψ1 = constructBC(u1, p)
    saveΨ(folderpsis, Ψ1, 0, p.params.Re)
    save_result(io, Ψ1, 0, 0, 0.0, p)
    u1 = [u1; p.params.Re / scl]

    p.params.Re = Re_start + ΔRe
    u2, _, _ = newton(f!, u0, p)
    Ψ2 = constructBC(u2, p)
    saveΨ(folderpsis, Ψ2, 1, p.params.Re)
    save_result(io, Ψ2, 1, 0, 0.0, p)
    u2 = [u2; p.params.Re / scl]

    # Step size norm
    s = norm(u2 - u1)

    for i in 2:(steps)
        time = @elapsed u, iter, _ = newton_for_continuation(f!, u1, u2, s, p)
        u, _, _ = newton_for_continuation(f!, u1, u2, s, p)

        u1 = u2
        u2 = u

        p.params.Re = u[end] * scl
        Ψ = constructBC(u[1:end-1], p)

        if isinteger(i / save_steps)
            saveΨ(folderpsis, Ψ, i, p.params.Re)
        end
        save_result(io, Ψ, i, iter, time, p)
    end
    close(io)
end

function saveΨ(foldername, Ψ, step, Re)
    writedlm("$(foldername)/psi_step$(@sprintf("%03d", step))_Re$(@sprintf("%07.3f", Re)).txt",  Ψ)
end

function save_result(io, Ψ, step, iter, time, p)
    @unpack Re, ic, i1, i2 = p.params 

    # Important: indices are transposed, mapping to physical space
    result = "$step,$Re,$(Ψ[ic,ic]),$(Ψ[i2,i2]),$(Ψ[i2,i1]),$(Ψ[i1,i2]),$(Ψ[i1,i1]),$iter,$(time)\n"
    write(io, result)
    flush(io)
end

function continuation_arclength_lsa(Ψi, p::CavityParameters, Re_start, Re_stop, ΔRe;
                                    verbose=false)
    @unpack n, scl, Ψ = p

    p.params.Re = Re_start
    @inbounds u0 = reshape(Ψi[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    u1, _, _ = newton(f!, u0, p)
    Ψ1 = constructBC(u1, p)
    save_result(io, Ψ1, 0, 0, 0.0, p)
    u1 = [u1; p.params.Re / scl]

    p.params.Re = Re_start + ΔRe
    u2, _, _ = newton(f!, u0, p)
    Ψ2 = constructBC(u2, p)
    u2 = [u2; p.params.Re / scl]

    # Step size norm
    s = norm(u2 - u1)


    if linearstability == true
        lambdas = Vector{typeof(u1[1:end-1])}(undef, steps)

    lsa_time1 = @elapsed lambdas[1] = linearstability_lambdas(u1[1:end-1], p)
    lsa_time2 = @elapsed lambdas[2] = linearstability_lambdas(u2[1:end-1], p)
    end

    if verbose == true
        if linearstability == true
            @printf("  %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n", "Step", "Re", "Iter.", "Ψc", "Newton", "λ1", "λ2", "λ3", "LSA")
            @printf("  %-6d %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f\n", 1, Re_series[1], 0, Ψ1[nxc, nyc], 0, lambdas[1][1], lambdas[1][1], lambdas[1][1], lsa_time1)
            @printf("  %-6d %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f\n", 2, Re_series[2], 0, Ψ2[nxc, nyc], 0, lambdas[2][1], lambdas[2][2], lambdas[2][3], lsa_time2)
        else
            @printf("  %-6s %-6s %-6s %-6s %-6s\n", "Step", "Re", "Iter.", "Ψc", "Newton")
            @printf("  %-6d %-6.2f %-6.2f %-6.2f %-6.2f\n", 1, Re_series[1], 0, Ψ1[nxc, nyc], 0)
            @printf("  %-6d %-6.2f %-6.2f %-6.2f %-6.2f\n", 2, Re_series[2], 0, Ψ2[nxc, nyc], 0)
        end
    end

    for i in 3:(steps)
        newton_time = @elapsed u, iter, _ = newton_for_continuation(f!, u1, u2, s, p)
        u, _, _ = newton_for_continuation(f!, u1, u2, s, p)

        u1 = u2
        u2 = u

        if linearstability == true
            lsa_time = @elapsed lambdas[i] = linearstability_lambdas(u[1:end-1], p)
        end

        Re_series[i] = u[end] * p.scl
        p.Re = Re_series[i]

        sol[i] = constructBC(u[1:(end - 1)], p)


        if verbose == true
            if linearstability == true
                @printf("  %-6d %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f\n", i, Re_series[i], iter, Ψ[nxc, nyc], newton_time, lambdas[i][1], lambdas[i][2], lambdas[i][3], lsa_time)
            else
                @printf("  %-6d %-6.2f %-6.2f %-6.2f %-6.2f\n", i, Re_series[i], iter, Ψ[nxc, nyc], newton_time)
            end
        end

    end

    if linearstability == true
        return sol, lambdas, Re_series
    end
    return sol, Re_series
end

function newton_for_continuation(f!::F, x1, x2, s, p::CavityStruct; tolmax = 1e-10,
                                 maxiter = 100) where {F}
    @unpack scl = p.params
    x = copy(x2)
    xi = x[1:(end - 1)]

    dim = size(x, 1) - 1

    eps = 1e-8

    v = x2 - x1
    xp = x2 + (s / norm(v)) * v

    fx = similar(x)
    fxi_Re = zeros(dim)

    J = zeros(dim, dim)
    J2 = zeros(dim + 1, dim + 1)

    dx = zeros(dim + 1)
    cache = FiniteDiff.JacobianCache(xi)

    iter = 0
    tol = 1.0
    while tol > tolmax && iter < maxiter
        p.params.Re = x[end] * scl # unscaled reynolds

        xi = x[1:(end - 1)]
        fxi = fx[1:(end - 1)]

        f!(fxi, xi, p)
        fx[1:(end - 1)] = fxi
        fx[end] = v' * (x - xp)

        FiniteDiff.finite_difference_jacobian!(J, (fxi, xi) -> f!(fxi, xi, p), xi, cache)
        # jacobian!(J, f!, xi, p; dx = eps)

        J2[1:dim, 1:dim] = J

        # Calculate change in "Reynolds variable"
        p.params.Re = (x[end] + eps) * scl # unscaled Reynolds
        f!(fxi_Re, xi, p)
        J2[1:dim, dim + 1] = (fxi_Re - fxi) / eps

        J2[dim + 1, :] = vec(v)

        dx = J2 \ (-fx)
        @. x = x + dx

        tol = norm(dx)
        iter += 1
    end

    return x, iter, tol
end
