function Cavity4Sided(mesh::ChebyshevMesh, reynolds::Number)
    # Neumann boundary conditions on the 2D sides default to zero
    bcleft = zeros(mesh.ny+1)
    bcright = zeros(mesh.ny+1)
    bcbottom = zeros(mesh.nx+1)
    bctop = zeros(mesh.nx+1)

    Cavity4Sided(mesh, bcleft, bcright, bcbottom, bctop, reynolds)
end

function constructBCmatrix(D::Matrix)
    n = size(D, 1) - 1
    M = [D[1,2] D[1,n]; D[n+1,2] D[n+1,n]]
    return inv(M)
end

function constructΨ(probl::Cavity4Sided, Ψi::Matrix)
    nx = probl.mesh.nx
    ny = probl.mesh.ny
    Dx1 = probl.mesh.diffx1
    Dy1 = probl.mesh.diffy1
    bcleft = probl.bcleft
    bcright = probl.bcright
    bcbottom = probl.bcbottom
    bctop = probl.bctop

    Ψ = zeros(nx+1, ny+1)
    Ψ[3:nx-1, 3:ny-1] = Ψi 

    # Setup helper matrices
    Mxinv = constructBCmatrix(Dx1) 
    Myinv = constructBCmatrix(Dy1) 

    h1 = Ψ[3:nx-1, 3:ny-1] * Dy1[1, 3:nx-1]
    h2 = Ψ[3:nx-1, 3:ny-1] * Dy1[ny+1, 3:ny-1]

    Ψ[3:nx-1, 2]  = Myinv[1, 1] * (bctop[3:nx-1] - h1) + Myinv[1, 2] * (bcbottom[3:nx-1] - h2)
    Ψ[3:nx-1, ny] = Myinv[2, 1] * (bctop[3:nx-1] - h1) + Myinv[2, 2] * (bcbottom[3:nx-1] - h2)

    h1 = Dx1[1, 3:nx-1]' * Ψ[3:nx-1, 2:ny]
    h2 = Dx1[nx+1, 3:nx-1]' * Ψ[3:nx-1, 2:ny]

    Ψ[2, 2:ny] = Mxinv[1, 1] * (bcright[2:ny]' - h1) + Mxinv[1, 2] * (bcleft[2:ny]' - h2)
    Ψ[nx, 2:ny] = Mxinv[2, 1] * (bcright[2:ny]' - h1) + Mxinv[2, 2] * (bcleft[2:ny]'- h2)

    return Ψ
end

function calculateΨinitial(probl::Cavity4Sided; initialguess="rand", Δt=1.::Number, nbtimesteps=200::Integer)
    nx = probl.mesh.nx
    ny = probl.mesh.ny

    if initialguess == "rand"
        Ψold = 1e-3*randn((nx+1),(nx+1))
    elseif initialguess == "zeros"
        Ψold = zeros((nx+1),(ny+1))
    else
        throw("initualguess $initialguess is not an option")
    end

    # Progress in time to get a stable solution
    for _ in 1:nbtimesteps
        # Implicit Euler to do one time step
        function ftimestep(ψi) 
            return rhstime(probl, Δt, Ψold, ψi)
        end
        ψiold = vec(Ψold[3:nx-1, 3:ny-1])
        ψi, _, _, _ = newton(ftimestep, ψiold)

        Ψi = reshape(ψi, (nx-3,nx-3))
        Ψold = CavityFlow.constructΨ(probl, Ψi)
    end
    return Ψold
end
