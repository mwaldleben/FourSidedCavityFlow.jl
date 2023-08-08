"""
    Ψsteady, iter, tol = steady_state(Ψ0, p; abstol = 1e-10, maxiters = 100)

Computes the steady-state solution for the four-sided cavity by employing the Newton algorithm.
Returns the calculated steady-state solution from a given initial guess `Ψ0` and 
the number of iterions `iter` and tolerance  `tol` of Newton's method.
"""
function steadystate(Ψ0, p::CavityStruct; abstol = 1e-10, maxiters = 100)
    @unpack n = p.params

    @inbounds u0 = reshape(Ψ0[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    Ψsteady = similar(Ψ0)

    u_steady, iter, tol = newton(f!, u0, p; abstol = abstol, maxiters = maxiters)

    @inbounds @views Ψsteady[3:(n - 1), 3:(n - 1)][:] .= u_steady
    Ψsteady = constructBC(u_steady, p)

    return Ψsteady, iter, tol
end
