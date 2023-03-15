@testset "solve.jl" begin
    @testset "jacobian" begin
        func(x) = @. x^2
        n = 3
        x = ones(n)
        Jref = 2 * I(n) 
        J = CavityFlow.jacobian(x, func)

        @test J ≈ Jref atol=1e-7
    end
    @testset "newton" begin
        func(x) = @. x^3
        n = 3
        x0 = ones(n)
        x, iter, tol, isconverged = CavityFlow.newton(func, x0)

        @test x ≈ zeros(n) atol=1e-7
    end
end
