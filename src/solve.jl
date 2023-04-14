function solve_steadystate(Ψ0, p::CavityParameters; tolmax = 1e-10, maxiter = 100)
    @unpack n = p

    @inbounds u0 = reshape(Ψ0[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    _, iter, tol = newton(f!, u0, p; tolmax = 1e-10, maxiter = 100)

    construct_BC!(p)

    return copy(p.Ψ), iter, tol
end

function solve_timestepping(Ψstart, p::CavityParameters, Δt, nb_timesteps)
    @unpack n, Ψ, Ψ0 = p

    @inbounds Ψ0 .= Ψstart
    @inbounds u0 = reshape(Ψstart[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    ft!(fu, u, p) = ftime!(fu, u, p, Δt)

    for step in 1:nb_timesteps
        u0, iter, tol = newton(ft!, u0, p)

        # println("Step $(step): $(iter) iterations")

        Ψ[3:(n - 1), 3:(n - 1)] .= reshape(u0, (n - 3, n - 3))

        Ψ0 .= Ψ
    end

    construct_BC!(p)

    return Ψ
end

function solve_timestepping_save(Ψstart, p::CavityParameters, Δt, steps)
    @unpack n, Ψ, Ψ0 = p

    @inbounds Ψ0 .= Ψstart
    @inbounds u0 = reshape(Ψstart[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    fu = similar(u0)

    ft!(fu, u, p) = ftime!(fu, u, p, Δt)

    sol = Vector{typeof(Ψstart)}(undef, steps)
    time_series = Vector{Float64}(undef, steps)

    sol[1] = Ψstart
    time_series[1] = 0

    for i in 1:(steps - 1)
        u0, _, _ = newton(ft!, u0, p)

        Ψ[3:(n - 1), 3:(n - 1)] .= reshape(u0, (n - 3, n - 3))
        construct_BC!(p)
        Ψ0 .= Ψ

        sol[i + 1] = Ψ
        time_series[i + 1] = Δt * i
    end

    return sol, time_series
end

function solve_continuation(Ψstart, p::CavityParameters, Re_start, ΔRe, steps)
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

    for i in 1:(steps - 2)
        u, _, _ = newton_continuation(f!, u1, u2, s, p)

        u1 = u2
        u2 = u

        Re_series[i + 2] = u[end] * p.scl
        p.Re = Re_series[i + 2]

        sol[i + 2] = construct_BC(p)
    end

    return sol, Re_series
end

function newton_continuation(f!::F, x1, x2, s, p; tolmax = 1e-10, maxiter = 100) where {F}
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

function newton(f!, x0, p; tolmax = 1e-10, maxiter = 100)
    x = copy(x0)
    dim = size(x, 1)

    fx = similar(x)
    f!(fx, x, p)

    J = zeros(dim, dim)
    dx = zeros(dim)
    cache = FiniteDiff.JacobianCache(x)

    iter = 0
    tol = 1.0

    f1!(fx, x) = f!(fx, x, p)

    while tol > tolmax && iter < maxiter
        FiniteDiff.finite_difference_jacobian!(J, f1!, x, cache)
        # jacobian!(J, f1!, x; dx = 1e-8)

        dx .= J \ (-fx)
        @. x = x + dx

        f!(fx, x, p)

        tol = norm(dx)
        iter += 1
    end

    return x, iter, tol
end

function jacobian!(J, f!, x; dx = 1e-8)
    dim = size(x, 1)

    Id = I(dim) * dx

    fx1 = similar(x)
    fx2 = similar(x)
    f!(fx1, x)

    @inbounds for m in 1:dim
        f!(fx2, x + Id[:, m])
        J[:, m] = @. (fx2 - fx1) / dx
    end
end

# Non allocating version
# function jacobian!(J, f!, x, xm, fx1, fx2; dx = 1e-8)
#     dim = size(x, 1)
#
#     f!(fx1, x)
#     xm .= x
#
#     @inbounds @simd for m in 1:dim
#         xm[m] = xm[m] + dx 
#
#         f!(fx2, xm)
#         J[:, m] = @. (fx2 - fx1) / dx 
#
#         xm[m] = xm[m] - dx 
#     end
# end
