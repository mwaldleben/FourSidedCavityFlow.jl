using NS2DBenchmarkSolver
using Test

@testset "diff.jl" begin
    n = 2
    Dₓ = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
    x = [-1; 0; 1]
    @test NS2DBenchmarkSolver.diffcheb(n)[1] ≈ x 
    @test NS2DBenchmarkSolver.diffcheb(n)[2] ≈ Dₓ 
end
