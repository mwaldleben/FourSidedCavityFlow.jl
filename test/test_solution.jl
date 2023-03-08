@testset "solution.jl" begin
    @testset "Solution1D" begin 
        x = [0; 1; 2]
        y = [1; 1; 1]
        span = (-1, 1)

        sol = Solution1D(x, y)

        @test sol.nodes == x
        @test sol.vals == y 
    end
    @testset "Solution1D plot recipes" begin
        @test skip=true
    end
end
