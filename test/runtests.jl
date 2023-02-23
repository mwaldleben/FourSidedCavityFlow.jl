using NS2DBenchmarkSolver
using Test

@testset "chebychev.jl" begin
    N = 0
    D = [0] 
    x = [1.] 
    @test NS2DBenchmarkSolver.chebychev(N) == (D, x)

    N = 2
    D = [1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5]
    x = [1.0; 0.0; -1.0]
    @test NS2DBenchmarkSolver.chebychev(N) .≈ (D, x)
end
