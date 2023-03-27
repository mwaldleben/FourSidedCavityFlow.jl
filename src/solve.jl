function solve_steadystate(Ψ0, p::CavityParameters)
    @unpack n = p

    @inbounds u0 = reshape(Ψ0[3:n-1,3:n-1],(n-3)*(n-3))

    _, iter, tol = newton(f!,u0,p)

    construct_BC!(p)

    return p.Ψ, iter, tol
end

function solve_timestepping(Ψstart, p::CavityParameters, Δt, nb_timesteps)
    @unpack n,Ψ,Ψ0 = p

    @inbounds Ψ0 .= Ψstart
    @inbounds u0 = reshape(Ψstart[3:n-1,3:n-1],(n-3)*(n-3))
    fu = similar(u0)

    ft!(fu,u,p) = ftime!(fu,u,p,Δt)

    for step in 1:nb_timesteps
        u0, iter, tol = newton(ft!,u0,p)

        Ψ[3:n-1,3:n-1] .= reshape(u0, (n-3,n-3))

        Ψ0 .= Ψ 

        time = [(Δt * step)]
    end

    construct_BC!(p)

    return Ψ
end

function newton(f!,x0, p; tolmax=1e-10, maxiter=100)
    x = copy(x0)
    dim = size(x, 1)

    fx = similar(x)
    f!(fx,x,p)

    J = zeros(dim,dim)
    dx = zeros(dim) 
    cache = FiniteDiff.JacobianCache(x)

    iter = 0
    tol = 1.0
    while tol>tolmax && iter<maxiter
        FiniteDiff.finite_difference_jacobian!(J,(fx,x) -> f!(fx,x,p),x,cache)
        # jacobian!(J,f!,x,p;dx=1e-8)

        dx .= J \ (-fx)
        @. x = x + dx        

        f!(fx,x,p)

        tol = norm(dx)
        iter += 1 
    end

    return x, iter, tol
end

function jacobian!(J,f!,x,p;dx=1e-8)
    dim = size(x, 1)

    Id = I(dim) * dx

    f1 = similar(x)
    f2 = similar(x)
    f!(f1,x,p)

    for m = 1:dim
        f!(f2,x+Id[:,m],p)
        J[:,m] = (f2 - f1) / Id[m,m]
    end
end
