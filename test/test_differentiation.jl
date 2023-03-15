@testset "differentation.jl" begin
    n = 2
    D1ref = [1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5]
    nodesref = [1; 0; -1]

    nodes, D1, _, _ = CavityFlow.diffchebyshev(n)
    @test nodes ≈ nodesref
    @test D1 ≈ D1ref

    # Test when using another reference length
    n = 2
    D1ref = [1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5]
    nodesref = [1; 0; -1]

    nodes, D1, _, _ = CavityFlow.diffchebyshev(n, length=2)
    @test nodes ≈ nodesref skip=true
    @test D1 ≈ D1ref skip=true
end
