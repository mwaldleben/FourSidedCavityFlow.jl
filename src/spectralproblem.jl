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

function ExamplePoisson1D(mesh::SpectralMesh1D, func::Function, (bcmin, bcmax))
    Dₓₓ  = mesh.diff2mat
    f = func(mesh.nodes)
    Dₓₓ, f = applyBC1D(Dₓₓ, f, bcmin, bcmax)
    return ExamplePoisson1D(mesh, Dₓₓ, f, bcmin, bcmax)
end

function solve(probl::ExamplePoisson1D)
    u = probl.diff2matBC \ probl.fBC

    if typeof(probl.bcmin) == BCDirichlet1D && typeof(probl.bcmax) == BCDirichlet1D 
        u = [0; u; 0]
    end

    sol = Solution1D(probl.mesh.nodes, u)
    return  sol
end
