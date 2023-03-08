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

function setBC2D(probl::Cavity4Sided, bcxmin::BCNeumann2D, bcxmax::BCNeumann2D, bcymin::BCNeumann2D, bcymax::BCNeumann2D)
    probl.bcxmin = bcxmin
    probl.bcxmax = bcxmax
    probl.bcymin = bcymin
    probl.bcymax = bcymax
end

function setBC2D(probl::Cavity4Sided, bcxmin::Function, bcxmax::Function, bcymin::Function, bcymax::Function)
    probl.bcxmin = BCNeumann2D(bcxmin(probl.mesh.ynodes))
    probl.bcxmax = BCNeumann2D(bcxmax(probl.mesh.ynodes))
    probl.bcymin = BCNeumann2D(bcymin(probl.mesh.xnodes))
    probl.bcymax = BCNeumann2D(bcymax(probl.mesh.xnodes))
end
