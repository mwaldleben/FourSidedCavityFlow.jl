@testset "continuation.jl" begin
    @testset "continuation_arclength" begin
        n = 6
        Re_start = 100
        p = CF.setup_struct(n, Re_start)
        Ψ0 = zeros((n + 1, n + 1))

        ΔRe = -1
        Re_steps = 5

        foldercont = "cont_test"
        mkdir(foldercont)
        sol, Re_series = CF.continuation_arclength(foldercont, Ψ0, p, Re_start, ΔRe, Re_steps)
        sol2, Re_series2 = load("$foldercont/psis.jld2", "sol", "Re_series")

        Ψ_ref = [0 0 0 0 0 0 0
            0 0.020364765599066 0.055490465905256 0.059511890031363 0.066044664922850 0.022395922147029 0
            0 -0.076815377070194 -0.086485356985073 -0.078484213580964 -0.031644130190074 -0.066261178052600 0
            0 -0.089424692036144 -0.150362055040675 -0.185477150143061 -0.150362055040676 -0.089424692036145 0
            0 -0.066261178052599 -0.031644130190072 -0.078484213580962 -0.086485356985071 -0.076815377070193 0
            0 0.022395922147029 0.066044664922850 0.059511890031364 0.055490465905256 0.020364765599066 0
            0 0 0 0 0 0 0]

        @test sol[end, :, :] ≈ Ψ_ref
        @test sol2[end, :, :] ≈ Ψ_ref

        rm(foldercont; recursive = true)
    end
    @testset "newton_continuation" begin
        n = 6
        Re = 100
        p = CF.setup_struct(n, Re)

        # Augmented system
        x1 = ones((n - 3) * (n - 3) + 1)
        x2 = 2 * ones((n - 3) * (n - 3) + 1)
        s = 0.05

        x, iter, tol = CF.newton_continuation(CF.f!, x1, x2, s, p; abstol = 1e-10,
            maxiters = 1)
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
end
