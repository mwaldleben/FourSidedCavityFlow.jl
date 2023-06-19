@testset "differentation.jl" begin
    n = 2
    D1ref = [1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5]
    nodesref = [1; 0; -1]

    nodes, D1, _, _ = CF.diff_chebyshev(n)
    @test nodes ≈ nodesref
    @test D1 ≈ D1ref
end
