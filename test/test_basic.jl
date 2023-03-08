@testset "differentation.jl" begin
    n = 2
    Dx1ref = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
    xref = [-1; 0; 1]

    x, Dx1, _, _ = NS2DBenchmarkSolver.diffchebychev(n)
    @test x ≈ xref
    @test Dx1 ≈ Dx1ref

    # Test if set to another intervall
    n = 2
    Dx1ref = [1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5]
    xref = [1; 0; -1]

    x, Dx1, _, _ = NS2DBenchmarkSolver.diffchebychev(n, span=(1, -1))
    @test x ≈ xref
    @test Dx1 ≈ Dx1ref
end

@testset "boundarycondition.jl" begin
    val = 1.1
    bc1 = BCDirichlet1D(val)
    bc2 = BCNeumann1D(val)
    @test bc1.val == val
    @test bc2.val == val
end

@testset "ExamplePoisson1D" begin
    nbcells = 16
    span = (-1, 1)

    mesh = SpectralMesh1D(nbcells, span)

    bcmin = BCDirichlet1D(0.0) 
    bcmax = BCDirichlet1D(0.0) 

    func(x) = exp(4*x) 
    probl = NS2DBenchmarkSolver.ExamplePoisson1D(mesh, func, (bcmin, bcmax)) 

    @test probl.bcmin.val == bcmin.val
    @test probl.bcmax.val == bcmax.val

    sol = NS2DBenchmarkSolver.solve(probl)

    exactsol(x) = @. (exp(4*x) - sinh(4)*x - cosh(4))/16
    exactvals = exactsol(mesh.nodes)

    @test sol.vals ≈ exactvals
end
