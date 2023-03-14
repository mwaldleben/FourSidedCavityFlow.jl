function solve(probl::Example1D)
    u = probl.diff2matBC \ probl.fBC

    if typeof(probl.bcmin) == BCDirichlet1D && typeof(probl.bcmax) == BCDirichlet1D 
        u = [0; u; 0]
    end

    sol = Solution1D(probl.mesh.nodes, u)
    return  sol
end

function solve(probl::Cavity4Sided, Ψinit::Matrix; tol::Real=1e-12, maxiter::Integer=100)
    nx = probl.mesh.nx
    ny = probl.mesh.ny
    Dx1 = probl.mesh.diffx1
    Dy1 = probl.mesh.diffy1
    bcxmin = probl.bcxmin
    bcxmax = probl.bcxmax
    bcymin = probl.bcymin
    bcymax = probl.bcymax

    ψinit = vec(Ψinit[3:nx-1, 3:ny-1])

    # Solve stationary equation using newton-raphson
    function fnewton(ψint) 
       return rhs(probl, ψint)
    end
    ψint, iter, tol, isconverged = newton(fnewton, ψinit; tol=tol, maxiter=maxiter)

    Ψint = reshape(ψint, (nx-3, ny-3))

    Ψ = constructΨboundary(Ψint, Dx1, Dy1, bcxmin, bcxmax, bcymin, bcymax)

    sol = Solution2D(probl.mesh.xnodes, probl.mesh.ynodes, Ψ, isconverged, tol, iter)

    return sol
end

function solve(probl::Cavity4Sided; tol::Real=1e-12, maxiter::Integer=100)
    Ψinit = zeros((probl.mesh.nx+1),(probl.mesh.ny+1))

    return solve(probl, Ψinit; tol=tol, maxiter=maxiter)
end

function rhs(probl::Cavity4Sided, ψint)
    nx = probl.mesh.nx
    ny = probl.mesh.ny
    Dx1 = probl.mesh.diffx1
    Dy1 = probl.mesh.diffy1
    Dx2 = probl.mesh.diffx2
    Dy2 = probl.mesh.diffy2
    Dx4 = probl.mesh.diffx4
    Dy4 = probl.mesh.diffy4
    bcxmin = probl.bcxmin
    bcxmax = probl.bcxmax
    bcymin = probl.bcymin
    bcymax = probl.bcymax

    Ψint = reshape(ψint, (nx-3, ny-3))
    Ψ = constructΨboundary(Ψint, Dx1, Dy1, bcxmin, bcxmax, bcymin, bcymax)

    biharmΨ = Dx4*Ψ +  Ψ*Dy4' + 2*(Dx2*Ψ)*Dy2'
    laplΨ = Dx2*Ψ + Ψ*Dy2'
    nonlinterm = (Dx1*Ψ).*(laplΨ*Dy1') - (Dx1*laplΨ).*(Ψ*Dy1')
    
    FΨ = (1/probl.reynolds)*biharmΨ - nonlinterm
    Fψint = vec(FΨ[3:nx-1, 3:ny-1])

    return Fψint
end

function rhstime(probl::Cavity4Sided, Δt::Real, Ψold::Matrix, ψint::Vector)
    nx = probl.mesh.nx
    ny = probl.mesh.ny
    Dx1 = probl.mesh.diffx1
    Dy1 = probl.mesh.diffy1
    Dx2 = probl.mesh.diffx2
    Dy2 = probl.mesh.diffy2
    Dx4 = probl.mesh.diffx4
    Dy4 = probl.mesh.diffy4
    bcxmin = probl.bcxmin
    bcxmax = probl.bcxmax
    bcymin = probl.bcymin
    bcymax = probl.bcymax

    Ψint = reshape(ψint, (nx-3, ny-3))
    Ψ = constructΨboundary(Ψint, Dx1, Dy1, bcxmin, bcxmax, bcymin, bcymax)

    Ψold = reshape(Ψold, (nx+1, ny+1))

    laplΨ = Dx2*Ψ + Ψ*Dy2'

    laplΨold = Dx2*Ψold + Ψold*Dy2'

    biharmΨ = Dx4*Ψ +  Ψ*Dy4' + 2*(Dx2*Ψ)*Dy2'

    nonlinterm = (Dx1*Ψ).*(laplΨ*Dy1') - (Dx1*laplΨ).*(Ψ*Dy1')
    
    FΨ = (1/probl.reynolds)*biharmΨ - nonlinterm

    FΨ = Δt*FΨ - laplΨ + laplΨold
    Fψint = vec(FΨ[3:nx-1, 3:ny-1])

    return Fψint
end

function jacobian(x, func::Function; dx=1e-8)
    dim = size(x,1)
    J = zeros((dim, dim))

    Id = I(dim) * dx

    Fx1 = func(x)
        
    for m = 1:dim
        Fx2 = func(x + Id[:, m])
        J[:, m] = (Fx2 - Fx1) / Id[m, m]
    end
    return J
end

function newton(func::Function, x0; tol=1e-12, maxiter=100)
    x = x0

    i = 0
    t = 1
    while t>tol && i<maxiter
        Fx = func(x)
        jac = jacobian(x, func)
        dx = jac \ (-Fx)
        x = x + dx        

        t = norm(dx)
        i += 1 
    end

    isconverged = true
    if i == maxiter 
        isconverged = false
        println("Newton method did not converge in $i iterations!")

    end

    return x, i, t, isconverged
end
