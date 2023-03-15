# Rectangular Mesh

abstract type RectangularMesh end

struct ChebyshevMesh <: RectangularMesh 
    nx::Int
    ny::Int
    lengthx::Real
    lengthy::Real
    xnodes::Vector
    ynodes::Vector
    diffx1::Matrix
    diffx2::Matrix
    diffx4::Matrix
    diffy1::Matrix
    diffy2::Matrix
    diffy4::Matrix
end


# Cavity problems

abstract type CavityProblem end

# Lid-driven 4 sided cavity flow problem 
# with Chebyshev discretization
mutable struct Cavity4Sided <: CavityProblem
    mesh::ChebyshevMesh
    bcleft::Vector
    bcright::Vector
    bcbottom::Vector
    bctop::Vector
    reynolds::Real
end
