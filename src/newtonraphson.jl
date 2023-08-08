"""
    x, iter, tol = newton1D(f::Function, x0, p; abstol = 1e-10, maxiters = 100)

1D Newton method. Returns the computed zero `x`, the number of iterations `iter` and the absolut tolerance `tol`.
p contains the the parameters for the function such that `fx = f(x, p)` 
"""
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
        @warn("Newton 1D did not converge in $iter iterations.")
    end

    return x, iter, tol
end

"""
    x, iter, tol = newton(f!, x0, p, newton_cache; abstol = 1e-10, maxiters = 100)

Newton-Raphson to compute the zero of a multiple-dimensional function `f`. An in-place version for the function
has to be provided `f(fx, x, p)`.

# Arguments
- `f!::Function`: In-place formulation of multi-dimensional unction, `f!(fx, x, p)`
- `x0::Real`: Inital guess
- `p`: Parameters and cache for the function
- `newton_cache::Real`: Tuple of variables used to cache for the Newton `newton_cache = (x, fx, dx, J, jac_cache)`, 
    `jac_cache` is created for the Jacobian calculation through the package `FiniteDiff.jl`
- `abstol::Real=1e-10`: Absolut tolerance  
- `maxiters::Integer=100`: Maximum iterations for the Newton algorithm 
"""
function newton(f!::Function, x0, p, newton_cache; abstol = 1e-10, maxiters = 100)
    x, fx, dx, J, jac_cache = newton_cache

    @inbounds x .= x0
    f!(fx, x, p)

    iter = 0
    tol = 1.0

    while tol > abstol && iter < maxiters
        FiniteDiff.finite_difference_jacobian!(J, (fx, x) -> f!(fx, x, p), x, jac_cache)

        dx = J \ fx
        @. x -= dx

        f!(fx, x, p)

        tol = norm(dx)
        iter += 1
    end

    if tol > abstol && iter == maxiters
        @warn("Newton did not converge in $iter iterations.")
    end

    return x, iter, tol
end

"""
    x, iter, tol = newton(f!, x0, p; abstol = 1e-10, maxiters = 100)

Newton-Raphson to compute the zero of a multiple-dimensional function `f`. An in-place version for the function
has to be provided `f(fx, x, p)`. This function version doesn't use a cache for the Newton and allocates at every call.

# Arguments
- `f!::Function`: In-place formulation of multi-dimensional function, `f!(fx, x, p)`
- `x0::Real`: Inital guess
- `p`: Parameters and cache for the function
- `abstol::Real=1e-10`: Absolut tolerance  
- `maxiters::Integer=100`: Maximum iterations for the Newton algorithm 
"""
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
