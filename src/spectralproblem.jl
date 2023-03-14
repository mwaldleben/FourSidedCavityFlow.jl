function Example1D(mesh::SpectralMesh1D, func::Function, (bcmin, bcmax))
    Dx2  = mesh.diff2
    f = func.(mesh.nodes)
    Dx2, f = applyBC1D(Dx2, f, bcmin, bcmax)
    return Example1D(mesh, Dx2, f, bcmin, bcmax)
end

function Cavity4Sided(mesh::SpectralMesh2D, reynolds)
    # Boundary conditions on the 2D sides default to zero
    bcxmin = BCNeumann2D(zeros(mesh.ny+1))
    bcxmax = BCNeumann2D(zeros(mesh.ny+1))
    bcymin = BCNeumann2D(zeros(mesh.nx+1))
    bcymax = BCNeumann2D(zeros(mesh.nx+1))

    Cavity4Sided(mesh, bcxmin, bcxmax, bcymin, bcymax, reynolds)
end

function constructBChelpermat(D)
    n = size(D, 1) - 1
    M = [D[1,2] D[1,n]; D[n+1,2] D[n+1,n]]
    return inv(M)
end

function constructΨboundary(Ψint, Dx1, Dy1, bcxmin::BCNeumann2D, bcxmax::BCNeumann2D, bcymin::BCNeumann2D, bcymax::BCNeumann2D)
    nx = size(Ψint, 1) + 3
    ny = size(Ψint, 2) + 3 

    Ψ = zeros(nx+1, ny+1)
    Ψ[3:nx-1, 3:ny-1] = Ψint 

    # Setup helper matrices
    Mxinv = constructBChelpermat(Dx1) 
    Myinv = constructBChelpermat(Dy1) 

    h1 = Ψ[3:nx-1, 3:ny-1] * Dy1[1, 3:nx-1]
    h2 = Ψ[3:nx-1, 3:ny-1] * Dy1[ny+1, 3:ny-1]

    Ψ[3:nx-1, 2]  = Myinv[1, 1] * (bcymin.vals[3:nx-1] - h1) + Myinv[1, 2] * (bcymax.vals[3:nx-1] - h2)
    Ψ[3:nx-1, ny] = Myinv[2, 1] * (bcymin.vals[3:nx-1] - h1) + Myinv[2, 2] * (bcymax.vals[3:nx-1] - h2)

    h1 = Dx1[1, 3:nx-1]' * Ψ[3:nx-1, 2:ny]
    h2 = Dx1[nx+1, 3:nx-1]' * Ψ[3:nx-1, 2:ny]

    Ψ[2, 2:ny] = Mxinv[1, 1] * (bcxmin.vals[2:ny]' - h1) + Mxinv[1, 2] * (bcxmax.vals[2:ny]' - h2)
    Ψ[nx, 2:ny] = Mxinv[2, 1] * (bcxmin.vals[2:ny]' - h1) + Mxinv[2, 2] * (bcxmax.vals[2:ny]'- h2)

    return Ψ
end

function calculateinitialguess(probl::Cavity4Sided; Δt=1, nbtimesteps=200)
    nx = probl.mesh.nx
    ny = probl.mesh.ny
    bcxmin = probl.bcxmin
    bcxmax = probl.bcxmax
    bcymin = probl.bcymin
    bcymax = probl.bcymax

    # First guess 
    # Ψ = zeros((probl.mesh.xnbcells+1),(probl.mesh.ynbcells+1))
    Ψold = 1e-3*randn((nx+1),(nx+1))

    # Progress in time to get closer to a stable solution
    for _ in 1:nbtimesteps
        # Function with implicit Euler to do one time step
        function ftimestep(ψint) 
            return rhstime(probl, Δt, Ψold, ψint)
        end
        ψintold = vec(Ψold[3:nx-1, 3:ny-1])
        ψint, _, _, _ = newton(ftimestep, ψintold)

        Ψint = reshape(ψint, (nx-3,nx-3))
        Ψold = NS2DBenchmarkSolver.constructΨboundary(Ψint, probl.mesh.diffx1, probl.mesh.diffy1, bcxmin, bcxmax, bcymin, bcymax)
    end
    return Ψold
end
