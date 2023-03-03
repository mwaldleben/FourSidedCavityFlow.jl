# Spectral Mesh

abstract type SpectralMesh end

struct SpectralMesh1D <: SpectralMesh
    nbcells::Int
    min::Real
    max::Real
    nodes::Vector
    diff1mat::Matrix
    diff2mat::Matrix
    diff4mat::Matrix
    spectralmethod::String
end

struct SpectralMesh2D <: SpectralMesh
    xnbcells::Int
    ynbcells::Int
    xmin::Real
    xmax::Real
    ymin::Real
    ymax::Real
    xnodes::Vector
    ynodes::Vector
    diffx1mat::Matrix
    diffx2mat::Matrix
    diffx4mat::Matrix
    diffy1mat::Matrix
    diffy2mat::Matrix
    diffy4mat::Matrix
    spectralmethod::String
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
mutable struct ExamplePoisson1D <: SpectralProblem
    mesh::SpectralMesh1D
    diff2matBC::Matrix
    fBC::Vector
    bcmin::BC1D
    bcmax::BC1D
end

# Lid-driven 4 sided cavity flow problem 
mutable struct Cavity4Sided <: SpectralProblem
    mesh::SpectralMesh2D
    psiinit::Matrix
end


# Solution 

abstract type Solution end

struct Solution1D <: Solution 
    nodes::Vector
    vals::Vector
end

