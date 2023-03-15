function solve(probl::Example1D)
    n = probl.mesh.n
    D2 = probl.mesh.diff2
    rhs = probl.rhs

    # Apply boundary conditions (Dirichlet),
    # is fixed for this example problem
    A = D2[2:n, 2:n]
    b = rhs[2:n]

    # Solve system
    u = A \ b 

    # Add zero boundary values to solution
    u = [0; u; 0]

    return  u
end

function solve(probl::Cavity4Sided, Ψinitial::Matrix; tolmax::Real=1e-12, maxiter::Integer=100)
    nx = probl.mesh.nx
    ny = probl.mesh.ny

    ψi = vec(Ψinitial[3:nx-1, 3:ny-1])

    # Solve stationary equation using newton-raphson
    function fnewton(ψi) 
       return rhs(probl, ψi)
    end
    ψi, iter, tol, isconverged = newton(fnewton, ψi; tolmax=tolmax, maxiter=maxiter)

    Ψi = reshape(ψi, (nx-3,ny-3))
    Ψ = constructΨ(probl, Ψi)

    return Ψ, iter, tol, isconverged
end

function solve(probl::Cavity4Sided; tolmax::Real=1e-12, maxiter::Integer=100)
    nx = probl.mesh.nx
    ny = probl.mesh.ny

    Ψi = zeros((nx+1),(ny+1))

    return solve(probl, Ψi; tolmax=tolmax, maxiter=maxiter)
end

function rhs(probl::Cavity4Sided, ψi)
    nx = probl.mesh.nx
    ny = probl.mesh.ny
    Dx1 = probl.mesh.diffx1
    Dy1 = probl.mesh.diffy1
    Dx2 = probl.mesh.diffx2
    Dy2 = probl.mesh.diffy2
    Dx4 = probl.mesh.diffx4
    Dy4 = probl.mesh.diffy4
    Re = probl.reynolds

    Ψi = reshape(ψi, (nx-3, ny-3))
    Ψ = constructΨ(probl, Ψi)

    biharmΨ = Dx4*Ψ +  Ψ*Dy4' + 2*(Dx2*Ψ)*Dy2'
    laplΨ = Dx2*Ψ + Ψ*Dy2'
    nonlinterm = (Dx1*Ψ).*(laplΨ*Dy1') - (Dx1*laplΨ).*(Ψ*Dy1')
    
    FΨ = (1/Re)*biharmΨ - nonlinterm
    Fψint = vec(FΨ[3:nx-1, 3:ny-1])

    return Fψint
end

function rhstime(probl::Cavity4Sided, Δt::Real, Ψold::Matrix, ψi::Vector)
    nx = probl.mesh.nx
    ny = probl.mesh.ny
    Dx1 = probl.mesh.diffx1
    Dy1 = probl.mesh.diffy1
    Dx2 = probl.mesh.diffx2
    Dy2 = probl.mesh.diffy2
    Dx4 = probl.mesh.diffx4
    Dy4 = probl.mesh.diffy4
    Re = probl.reynolds

    Ψi = reshape(ψi, (nx-3, ny-3))
    Ψ = constructΨ(probl, Ψi)

    Ψold = reshape(Ψold, (nx+1, ny+1))

    laplΨ = Dx2*Ψ + Ψ*Dy2'
    laplΨold = Dx2*Ψold + Ψold*Dy2'

    biharmΨ = Dx4*Ψ +  Ψ*Dy4' + 2*(Dx2*Ψ)*Dy2'
    nonlinterm = (Dx1*Ψ).*(laplΨ*Dy1') - (Dx1*laplΨ).*(Ψ*Dy1')
    
    FΨ = (1/Re)*biharmΨ - nonlinterm

    FΨ = Δt*FΨ - laplΨ + laplΨold
    Fψint = vec(FΨ[3:nx-1, 3:ny-1])

    return Fψint
end

function jacobian(x::Vector, func::Function; dx=1e-8::Number)
    dim = size(x, 1)
    J = zeros((dim, dim))

    Id = I(dim) * dx

    Fx1 = func(x)
        
    for m = 1:dim
        Fx2 = func(x + Id[:, m])
        J[:, m] = (Fx2 - Fx1) / Id[m, m]
    end
    return J
end

function newton(func::Function, x0::Vector; tolmax=1e-12::Number, maxiter=100::Integer)
    x = x0

    iter = 0
    tol = 1
    while tol>tolmax && iter<maxiter
        Fx = func(x)
        jac = jacobian(x, func)
        dx = jac \ (-Fx)
        x = x + dx        

        tol = norm(dx)
        iter += 1 
    end

    isconverged = true
    if iter == maxiter 
        isconverged = false
        println("Newton method did not converge in $iter iterations!")

    end

    return x, iter, tol, isconverged
end
