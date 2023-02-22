using NS2DBenchmarkSolver
using Test

@testset "NS2DBenchmarkSolver.jl" begin
    x = 2
    y = 2
    @test NS2DBenchmarkSolver.sum_values(x,y) == 4
end
