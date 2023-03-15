@testset "Example1D" begin
    n = 16
    length = 1

    mesh = SpectralMesh1D(n)

    func(x) = exp(4*x) 
    probl = Example1D(mesh, func) 

    sol = NS2DBenchmarkSolver.solve(probl)

    fexact(x) = @. (exp(4*x) - sinh(4)*x - cosh(4))/16
    solexact = fexact(mesh.nodes)

    @test sol.vals ≈ solexact 
end
