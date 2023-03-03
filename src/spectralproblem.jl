function ExamplePoisson1D(mesh::SpectralMesh1D, func::Function, (bcmin, bcmax))
    Dx2  = mesh.diff2mat
    f = func.(mesh.nodes)
    Dx2, f = applyBC1D(Dx2, f, bcmin, bcmax)
    return ExamplePoisson1D(mesh, Dx2, f, bcmin, bcmax)
end

function solve(probl::ExamplePoisson1D)
    u = probl.diff2matBC \ probl.fBC

    if typeof(probl.bcmin) == BCDirichlet1D && typeof(probl.bcmax) == BCDirichlet1D 
        u = [0; u; 0]
    end

    sol = Solution1D(probl.mesh.nodes, u)
    return  sol
end

function Cavity4Sided(mesh::SpectralMesh2D)
    nx = mesh.xnbcells
    ny = mesh.ynbcells
    psiinit = zeros(nx+1, ny+1)
    Cavity4Sided(mesh, psiinit)
end

function constructBChelpermat(D)
    n = size(D, 1) - 1
    M = [D[1,2] D[1,n]; D[n+1,2] D[n+1,n]]
    return inv(M)
end

function fillΨint(nx, ny; initialΨval=0)
    Ψ = zeros(nx+1,nx+1);

    if initialΨval ≢ 0
        Ψ[3:nx-1,3:ny-1] = fill(initialΨval, (nx-3, ny-3))
    end
    return  Ψ 
end

function constructΨboundary(Ψ, Dx1, Dy1, bcxmin::BCNeumann2D, bcxmax::BCNeumann2D, bcymin::BCNeumann2D, bcymax::BCNeumann2D)
    nx = size(Ψ, 1) - 1
    ny = size(Ψ, 2) - 1 

    # Setup helper matrices
    Mxinv = constructBChelpermat(Dx1) 
    Myinv = constructBChelpermat(Dy1) 

    h1 = Ψ[3:nx-1, 3:ny-1] * Dy1[1, 3:nx-1]'
    h2 = Ψ[3:nx-1, 3:ny-1] * Dy1[ny+1, 3:ny-1]'

    Ψ[3:nx-1, 2]  = Mxinv[1, 1] * (bcymin.vals[3:nx-1] - h1) + Myinv[1, 2] * (bcymax.vals[3:nx-1] - h2)
    Ψ[3:nx-1, ny] = Myinv[2, 1] * (bcymin.vals[3:nx-1] - h1) + Myinv[2, 2] * (bcymax.vals[3:nx-1] - h2)

    h1 = Dx1[1, 3:nx-1] * Ψ[3:nx-1, 2:ny]
    h2 = Dx1[nx+1, 3:nx-1] * Ψ[3:nx-1, 2:ny]

    Ψ[2, 2:ny] = Mxinv[1, 1] * (bcxmin.vals[2:ny]' - h1) + Mxinv[1, 2] * (bcxmax.vals[2:ny]' - h2)
    Ψ[nx, 2:ny] = Mxinv[2, 1] * (bcxmin.vals[2:ny]' - h1) + Mxinv[2, 2] * (bcxmax.vals[2:ny]'- h2)

    return Ψ
end

function solve(probl::Cavity4Sided)
    throw("unimplemented")
    return
end
