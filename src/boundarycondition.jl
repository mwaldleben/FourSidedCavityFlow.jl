function applyBC1D(A::Matrix, b::Vector, bcmin::BC1D, bcmax::BC1D)
    n = length(b)

    if typeof(bcmin) == BCDirichlet1D && typeof(bcmax) == BCDirichlet1D 
        A = A[2:n-1, 2:n-1] # pure Dirichlet boundary conditions 
        b = b[2:n-1]
        return A, b
    else
        throw("unimplemented")
    end
end

function BCNeumann2D(nodes, bcfunc::Function)
    vals = bcfunc.(nodes)
    return BCNeumann2D(vals)
end

function applyBC2D(probl::Cavity4Sided, bcxmin::BCNeumann2D, bcxmax::BCNeumann2D, bcymin::BCNeumann2D, bcymax::BCNeumann2D; initialΨval=0)
    nx = probl.mesh.xnbcells
    ny = probl.mesh.ynbcells 
    Dx1 = probl.mesh.diffx1mat
    Dy1 = probl.mesh.diffy1mat

    # Fill interior from initial value of Ψ
    Ψ = fillΨint(nx, ny; initialΨval)
    
    # Construct Ψ at the boundary 
    Ψ = constructΨboundary(Ψ, Dx1, Dy1, bcxmin, bcxmax, bcymin, bcymax)

    probl.psiinit = Ψ
end


function setBC2D(probl::Cavity4Sided, bcxmin::BCNeumann2D, bcxmax::BCNeumann2D, bcymin::BCNeumann2D, bcymax::BCNeumann2D; initialΨval=0)
    probl.bcxmin = bcxmin
    probl.bcxmax = bcxmax
    probl.bcymin = bcymin
    probl.bcymax = bcymax
end
