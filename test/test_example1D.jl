@testset "boundarycondition.jl" begin
    val = 1.1
    bc1 = BCDirichlet1D(val)
    bc2 = BCNeumann1D(val)
    @test bc1.val == val
    @test bc2.val == val
end

@testset "Example1D" begin
    n = 16
    length = 1

    mesh = SpectralMesh1D(n, length)

    bcmin = BCDirichlet1D(0.0) 
    bcmax = BCDirichlet1D(0.0) 

    func(x) = exp(4*x) 
    probl = Example1D(mesh, func, (bcmin, bcmax)) 

    @test probl.bcmin.val == bcmin.val
    @test probl.bcmax.val == bcmax.val

    sol = NS2DBenchmarkSolver.solve(probl)

    exactsol(x) = @. (exp(4*x) - sinh(4)*x - cosh(4))/16
    exactvals = exactsol(mesh.nodes)

    @test sol.vals ≈ exactvals
end
