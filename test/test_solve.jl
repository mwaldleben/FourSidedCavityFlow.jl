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
        steps = 7
        sol, Re_series = CavityFlow.solve_continuation(Ψ0, p, Re_start, ΔRe, steps)

        Ψ_ref = [0 0 0 0 0 0 0
                 0 0.022395922147029 0.066044664922850 0.059511890031364 0.055490465905256 0.020364765599066 0
                 0 -0.066261178052599 -0.031644130190072 -0.078484213580962 -0.086485356985071 -0.076815377070193 0
                 0 -0.089424692036144 -0.150362055040675 -0.185477150143061 -0.150362055040676 -0.089424692036145 0
                 0 -0.076815377070194 -0.086485356985073 -0.078484213580964 -0.031644130190074 -0.066261178052600 0
                 0 0.020364765599066 0.055490465905256 0.059511890031363 0.066044664922850 0.022395922147029 0
                 0 0 0 0 0 0 0]

        @test sol[steps] ≈ Ψ_ref
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
    @testset "newton_linearstability" begin
        n = 8
        dim = (n - 3) * (n - 3)
        Re0 = 66
        p = CavityFlow.setup_params(n, Re0)

        Ψ0 = zeros(n + 1, n + 1)
        Ψsteady, iter, tol = CavityFlow.solve_steadystate(Ψ0, p)

        u0 = Ψsteady[3:(n - 1), 3:(n - 1)][:]

        Re, iter, tol = CavityFlow.newton1D_linearstability(Re0, u0, p; tolmax = 1e-8,
                                                            maxiter = 20)

        Re_ref = 68.943470270431973

        @test Re≈Re_ref atol=1e-8
    end
    @testset "newton1D" begin
        function f(x, p)
            return p * x^3
        end

        p = 1
        x0 = 1.0
        x, iter, tol = CavityFlow.newton1D(f, x0, p)

        @test x≈0.0 atol=1e-7
    end
    @testset "newton" begin
        function f!(fx, x, p)
            @. fx = p * x^3
        end

        n = 3
        p = 1
        x0 = ones(n)
        x, iter, tol = CavityFlow.newton(f!, x0, p)

        @test x≈zeros(n) atol=1e-7
    end
    @testset "jacobian!" begin
        function f!(fx, x)
            @. fx = x^2
        end

        n = 3
        p = 1
        x = ones(n)
        J_ref = 2 * I(n)
        J = zeros(n, n)
        CavityFlow.jacobian!(J, f!, x)

        @test J≈J_ref atol=1e-7
    end
end
