@testset "rectangularmesh.jl" begin
    @testset "ChebyshevMesh" begin
        n = (2, 2)
        D1ref = [1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5]
        nodesref = [1; 0; -1]

        mesh = ChebyshevMesh(n)

        @test mesh.xnodes ≈ nodesref
        @test mesh.ynodes ≈ nodesref
        @test mesh.diffx1 ≈ D1ref
        @test mesh.diffy1 ≈ D1ref
    end
end
