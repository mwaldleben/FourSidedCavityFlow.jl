@testset "newtonraphson.jl" begin
    @testset "newton1D" begin
        function f(x, p)
            return p * x^3
        end

        p = 1
        x0 = 1.0
        x, iter, tol = CF.newton1D(f, x0, p)

        @test x≈0.0 atol=1e-7
    end
    @testset "newton" begin
        function f!(fx, x, p)
            @. fx = p * x^3
        end

        n = 3
        p = 1
        x0 = ones(n)

        x, iter, tol = CF.newton(f!, x0, p)

        @test x≈zeros(n) atol=1e-7
    end
end
