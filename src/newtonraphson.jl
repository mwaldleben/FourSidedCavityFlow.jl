function newton1D(f, x0, p, tolmax = 1e-10, maxiter = 100)
    x = x0

    iter = 0
    tol = 1.0

    while tol > tolmax && iter < maxiter
        fx = f(x, p)

        # Fixed step size!
        x1 = x + 1e-8
        fx1 = f(x1, p)
        dfx = (fx1 - fx) / 1e-8

        dx = -fx / dfx
        x = x + dx

        tol = abs(dx)
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
