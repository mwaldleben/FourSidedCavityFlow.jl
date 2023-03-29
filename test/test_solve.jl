@testset "solve.jl" begin
    @testset "solve_steadystate" begin
        n = 6
        Re = 100
        p = CavityFlow.setup_params(n, Re)

        Ψ0 = zeros(n + 1, n + 1)
        Ψsteady, iter, tol = CavityFlow.solve_steadystate(Ψ0, p)

        Ψsteady_ref = [0 0 0 0 0 0 0
                       0 0.0223534291102081 0.0655941971343996 0.0589510766581025 0.0555732229353761 0.0204248917264253 0
                       0 -0.0664329942047326 -0.0337131774973950 -0.0810657532138882 -0.0857836868595323 -0.0764539684037562 0
                       0 -0.0892424343792270 -0.150714381583066 -0.187874748844775 -0.150714381583067 -0.0892424343792272 0
                       0 -0.0764539684037564 -0.0857836868595334 -0.0810657532138900 -0.0337131774973966 -0.0664329942047329 0
                       0 0.0204248917264252 0.0555732229353758 0.0589510766581022 0.0655941971343994 0.0223534291102081 0
                       0 0 0 0 0 0 0]

        @test Ψsteady ≈ Ψsteady_ref
    end
    @testset "solve_timestepping" begin
        n = 6
        Re = 100
        p = CavityFlow.setup_params(n, Re)

        Δt = 1
        steps = 200
        Ψ0 = ones((n + 1, n + 1))

        Ψ200 = CavityFlow.solve_timestepping(Ψ0, p, Δt, steps)

        Ψ200_ref = [0 0 0 0 0 0 0
                    0 0.0223534291102081 0.0655941971343996 0.0589510766581025 0.0555732229353761 0.0204248917264253 0
                    0 -0.0664329942047326 -0.0337131774973950 -0.0810657532138882 -0.0857836868595323 -0.0764539684037562 0
                    0 -0.0892424343792270 -0.150714381583066 -0.187874748844775 -0.150714381583067 -0.0892424343792272 0
                    0 -0.0764539684037564 -0.0857836868595334 -0.0810657532138900 -0.0337131774973966 -0.0664329942047329 0
                    0 0.0204248917264252 0.0555732229353758 0.0589510766581022 0.0655941971343994 0.0223534291102081 0
                    0 0 0 0 0 0 0]

        @test Ψ200 ≈ Ψ200_ref
    end
    @testset "solve_continuation" begin
        n = 6
        Re = 100
        p = CavityFlow.setup_params(n, Re)

        Ψ0 = zeros((n + 1, n + 1))

        Re_start = Re
        ΔRe = -1
        steps = 2
        CavityFlow.solve_continuation(Ψ0, p, Re_start, ΔRe, steps)
    end
    @testset "newton_continuation" begin
        n = 6
        Re = 100
        p = CavityFlow.setup_params(n, Re)

        # augmented system
        x1 = ones((n - 3) * (n - 3) + 1)
        x2 = 2 * ones((n - 3) * (n - 3) + 1)
        s = 0.05

        x, iter, tol = CavityFlow.newton_continuation(CavityFlow.f!, x1, x2, s, p;
                                                      tolmax = 1e-10, maxiter = 1)

        x_ref = [
            1.81433267406029,
            1.82898296631620,
            1.81433272036535,
            1.80516123907409,
            1.81997644436709,
            1.80516123057598,
            1.81433271763541,
            1.82898296504399,
            1.81433267684793,
            3.81251824872209,
        ]
        @test x≈x_ref atol=1e-2
    end
    @testset "newton" begin
        function fnewton!(fx, x, p)
            @. fx = p * x^3
        end

        n = 3
        p = 1
        x0 = ones(n)
        x, iter, tol = CavityFlow.newton(fnewton!, x0, p)

        @test x≈zeros(n) atol=1e-7
    end
    @testset "jacobian!" begin
        function fjac!(fx, x, p)
            @. fx = p * x^2
        end

        n = 3
        p = 1
        x = ones(n)
        J_ref = 2 * I(n)
        J = zeros(n, n)
        CavityFlow.jacobian!(J, fjac!, x, p)

        @test J≈J_ref atol=1e-7
    end
end
