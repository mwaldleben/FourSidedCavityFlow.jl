@testset "spectralmesh.jl" begin
    @testset "SpectralMesh1D" begin
        nbcells = 2
        span = (-1, 1)
        Dx1ref = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
        xref = [-1; 0; 1]

        mesh = SpectralMesh1D(nbcells, span)

        @test mesh.nodes ≈ xref
        @test mesh.diff1mat ≈ Dx1ref
        @test mesh.spectralmethod == "chebychev"
    end
    @testset "SpectralMesh2D" begin
        nbcells = (2, 2)
        xspan = (-1, 1)
        yspan = (-1, 1)
        Dx1ref = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
        xref = [-1; 0; 1]

        mesh = SpectralMesh2D(nbcells, xspan, yspan)

        @test mesh.xnodes ≈ xref
        @test mesh.ynodes ≈ xref
        @test mesh.diffx1mat ≈ Dx1ref
        @test mesh.diffy1mat ≈ Dx1ref
        @test mesh.spectralmethod == "chebychev"
    end
end
