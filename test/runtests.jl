using NS2DBenchmarkSolver
using Test

@testset "differentation.jl" begin
    n = 2
    Dₓ = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
    x = [-1; 0; 1]
    @test NS2DBenchmarkSolver.diffchebychev(n)[1] ≈ x
    @test NS2DBenchmarkSolver.diffchebychev(n)[2] ≈ Dₓ
end

@testset "boundarycondition.jl" begin
    val = 1.1
    bc1 = NS2DBenchmarkSolver.BCDirichlet1D(val)
    bc2 = NS2DBenchmarkSolver.BCNeumann1D(val)
    @test bc1.val == val
    @test bc2.val == val
end

@testset "solution.jl" begin
    @testset "Solution1D" begin 
        x = [0; 1; 2]
        y = [1; 1; 1]

        sol = NS2DBenchmarkSolver.Solution1D(x, y)

        @test sol.nodes == x
        @test sol.vals == y 
    end
    @testset "Solution1D plot recipes" begin
        @test skip=true
    end
end

@testset "spectralmesh.jl" begin
    nbcells = 2
    coordmin = -1
    coordmax = 1
    Dₓ = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
    x = [-1; 0; 1]

    mesh = NS2DBenchmarkSolver.SpectralMesh1D(nbcells, (coordmin, coordmax))

    @test mesh.nodes ≈ x
    @test mesh.diff1mat ≈ Dₓ
    @test mesh.spectralmethod == "chebychev"
end

@testset "spectralproblem.jl" begin
    @testset "ExamplePoisson1D" begin
        nbcells = 16
        coordmin = -1
        coordmax = 1

        mesh = NS2DBenchmarkSolver.SpectralMesh1D(nbcells, (coordmin, coordmax))

        bcmin = NS2DBenchmarkSolver.BCDirichlet1D(0.0) 
        bcmax = NS2DBenchmarkSolver.BCDirichlet1D(0.0) 

        func(x) = exp.(4*x) 
        probl = NS2DBenchmarkSolver.ExamplePoisson1D(mesh, func, (bcmin, bcmax)) 

        @test probl.bcmin.val == bcmin.val
        @test probl.bcmax.val == bcmax.val

        sol = NS2DBenchmarkSolver.solve(probl)

        exactsol(x) = @. (exp(4*x) - sinh(4)*x - cosh(4))/16
        exactvals = exactsol(mesh.nodes)

        @test sol.vals ≈ exactvals
    end
end
