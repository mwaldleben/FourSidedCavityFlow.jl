# Spectral Mesh

abstract type SpectralMesh end

struct SpectralMesh1D <: SpectralMesh
    n::Int
    length::Real
    nodes::Vector
    diff1::Matrix
    diff2::Matrix
    diff4::Matrix
end

struct SpectralMesh2D <: SpectralMesh
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


# Boundary condition

abstract type BoundaryCondition end

abstract type BC1D <: BoundaryCondition end

struct BCDirichlet1D <: BC1D 
    val::Real
end

struct BCNeumann1D <: BC1D 
    val::Real
end

abstract type BC2D <: BoundaryCondition end

struct BCDirichlet2D <: BC2D
    vals::Vector
end

struct BCNeumann2D <: BC2D
    vals::Vector
end


# Spectral problems

abstract type SpectralProblem end

# Solve linear 1D BVP with Dirichlet or Neumann boundary conditions
# Example: u_xx = f(x) 
# Corresponding Matlab code, Program 13 and 33 in Trefethen
mutable struct Example1D <: SpectralProblem
    mesh::SpectralMesh1D
    diff2matBC::Matrix
    fBC::Vector
    bcmin::BC1D
    bcmax::BC1D
end

# Lid-driven 4 sided cavity flow problem 
mutable struct Cavity4Sided <: SpectralProblem
    mesh::SpectralMesh2D
    bcxmin::BCNeumann2D
    bcxmax::BCNeumann2D
    bcymin::BCNeumann2D
    bcymax::BCNeumann2D
    reynolds::Real
end


# Solution 

abstract type Solution end

struct Solution1D <: Solution
    nodes::Vector
    vals::Vector
end

struct Solution2D <: Solution
    xnodes::Vector
    ynodes::Vector
    vals::Matrix
    isconverged::Bool
    tol::Real
    iter::Integer
end

