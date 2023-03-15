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


# Spectral problems

abstract type SpectralProblem end

# Example problem: u_xx = f(x), u(-1) = u(1) = 0
# Solve linear 1D Boundary value problem with Dirichlet boundary conditions
# Corresponding Matlab code of program 13 in
# Spectral Methods in Matlab, Lloyd N. Trefethen
mutable struct Example1D <: SpectralProblem
    mesh::SpectralMesh1D
    rhs::Vector
end

# Lid-driven 4 sided cavity flow problem 
mutable struct Cavity4Sided <: SpectralProblem
    mesh::SpectralMesh2D
    bcleft::Vector
    bcright::Vector
    bcbottom::Vector
    bctop::Vector
    reynolds::Real
end
