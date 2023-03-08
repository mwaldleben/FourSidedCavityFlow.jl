function solve(probl::Cavity4Sided; tol::Real=1e-8, maxiter::Integer=100)
    nx = probl.mesh.xnbcells
    ny = probl.mesh.ynbcells
    Dx1 = probl.mesh.diffx1mat
    Dy1 = probl.mesh.diffy1mat
    bcxmin = probl.bcxmin
    bcxmax = probl.bcxmax
    bcymin = probl.bcymin
    bcymax = probl.bcymax

    ψinit = zeros((nx-3)*(ny-3))

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

function rhs(probl::Cavity4Sided, psiint)
    nx = probl.mesh.xnbcells
    ny = probl.mesh.ynbcells
    Dx1 = probl.mesh.diffx1mat
    Dy1 = probl.mesh.diffy1mat
    Dx2 = probl.mesh.diffx2mat
    Dy2 = probl.mesh.diffy2mat
    Dx4 = probl.mesh.diffx4mat
    Dy4 = probl.mesh.diffy4mat
    bcxmin = probl.bcxmin
    bcxmax = probl.bcxmax
    bcymin = probl.bcymin
    bcymax = probl.bcymax

    Ψint = reshape(psiint, (nx-3, ny-3))
    Ψ = constructΨboundary(Ψint, Dx1, Dy1, bcxmin, bcxmax, bcymin, bcymax)

    biharmΨ = Dx4*Ψ +  Ψ*Dy4' + 2*(Dx2*Ψ)*Dy2'
    laplΨ = Dx2*Ψ + Ψ*Dy2'
    nonlinterm = (Dx1*Ψ).*(laplΨ*Dy1') - (Dx1*laplΨ).*(Ψ*Dy1')
    
    FΨ = (1/probl.reynolds)*biharmΨ - nonlinterm
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

function newton(func::Function, x0; tol=1e-8, maxiter=100)
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
