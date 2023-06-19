function newton1D(f::Function, x0, p, abstol = 1e-10, maxiters = 100)
    x = x0

    iter = 0
    tol = 1.0

    while tol > abstol && iter < maxiters
        fx = f(x, p)

        # Fixed step size!
        x1 = x + 1e-8
        fx1 = f(x1, p)
        dfx = (fx1 - fx) / 1e-8

        dx = fx / dfx
        x -= dx

        tol = abs(dx)
        iter += 1
    end

    if tol > abstol && iter == maxiters
        @warn "Newton did not converge in $iter iterations."
    end

    return x, iter, tol
end

function newton(f!::Function, x0, p, newton_cache; abstol = 1e-10, maxiters = 100)
    x, fx, dx, J, jac_cache = newton_cache

    f1!(fx, x) = f!(fx, x, p)

    @inbounds x .= x0
    f1!(fx, x)

    iter = 0
    tol = 1.0

    while tol > abstol && iter < maxiters
        FiniteDiff.finite_difference_jacobian!(J, f1!, x, jac_cache)

        dx = J \ fx
        @. x -= dx

        f!(fx, x, p)

        tol = norm(dx)
        iter += 1
    end

    if tol > abstol && iter == maxiters
        @warn "Newton did not converge in $iter iterations."
    end

    return x, iter, tol
end

function newton(f!::Function, x0, p; abstol = 1e-10, maxiters = 100)
    x = similar(x0)
    fx = similar(x0)
    dx = similar(x0)

    dim = size(x0, 1)
    J = zeros(dim, dim)
    jac_cache = FiniteDiff.JacobianCache(x0)

    newton_cache = (x, fx, dx, J, jac_cache)
    x, iter, tol = newton(f!, x0, p, newton_cache; abstol = abstol, maxiters = maxiters)

    return x, iter, tol
end
