abstract type BoundaryCondition end

abstract type BC1D <: BoundaryCondition end

struct BCDirichlet1D <: BC1D 
    val::Real
end

struct BCNeumann1D <: BC1D 
    val::Real
end

function applyBC1D(A::Matrix, b::Vector, bcmin::BC1D, bcmax::BC1D)
    n = length(b)

    if typeof(bcmin) == BCDirichlet1D && typeof(bcmax) == BCDirichlet1D 
        A = A[2:n-1,2:n-1] # pure Dirichlet boundary conditions 
        b = b[2:n-1]
        return A, b
    else
        throw("unimplemented")
    end
end
