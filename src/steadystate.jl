function steadystate(Ψ0, p::CavityParameters; tolmax = 1e-10, maxiter = 100)
    @unpack n = p

    @inbounds u0 = reshape(Ψ0[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    _, iter, tol = newton(f!, u0, p; tolmax = 1e-10, maxiter = 100)

    construct_BC!(p)

    return copy(p.Ψ), iter, tol
end
