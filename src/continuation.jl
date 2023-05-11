function continuation_arclength(Ψstart, p::CavityParameters, Re_start, ΔRe, steps;
                                linearstability = false, verbose=false)
    @unpack n, scl, Ψ = p

    p.Re = Re_start
    @inbounds Ψ .= Ψstart
    @inbounds u0 = reshape(Ψ[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    u1, _, _ = newton(f!, u0, p)
    u1 = [u1; p.Re / scl]
    Ψ1 = construct_BC(p)

    p.Re = Re_start + ΔRe
    u2, _, _ = newton(f!, u0, p)
    u2 = [u2; p.Re / scl]
    Ψ2 = construct_BC(p)

    s = norm(u2 - u1)

    sol = Vector{typeof(Ψstart)}(undef, steps)
    Re_series = Vector{Float64}(undef, steps)

    sol[1] = Ψ1
    sol[2] = Ψ2

    Re_series[1] = Re_start
    Re_series[2] = Re_start + ΔRe

    nxc = Int(floor(n / 2) + 1)
    nyc = Int(floor(n / 2) + 1)

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

        @inbounds @views p.Ψ[3:(n - 1), 3:(n - 1)][:] .= u[1:(end - 1)]
        sol[i] = construct_BC(p)


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

function newton_for_continuation(f!::F, x1, x2, s, p; tolmax = 1e-10,
                                 maxiter = 100) where {F}
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
        p.Re = x[end] * p.scl # unscaled reynolds

        xi = x[1:(end - 1)]
        fxi = fx[1:(end - 1)]

        f!(fxi, xi, p)
        fx[1:(end - 1)] = fxi
        fx[end] = v' * (x - xp)

        FiniteDiff.finite_difference_jacobian!(J, (fxi, xi) -> f!(fxi, xi, p), xi, cache)
        # jacobian!(J, f!, xi, p; dx = eps)

        J2[1:dim, 1:dim] = J

        # Calculate change in "Reynolds variable"
        p.Re = (x[end] + eps) * p.scl # unscaled Reynolds
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
