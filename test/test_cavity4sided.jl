@testset "Cavity4Sided" begin
    @testset "Cavity4Sided methods" begin
        # Test constructBCMatrix function
        #nodes = [1; 1/2; -1/2; -1] 
        D = [19/6 -4 4/3 -1/2; 1 -1/3 -1 1/3; -1/3 1 1/3 -1; 1/2 -4/3 4 -19/6]

        Minvref = [-9 3; -3 9]/32
        Minv = NS2DBenchmarkSolver.constructBCmatrix(D)
        @test Minv ≈ Minvref

        # Test boundary reconstruction Ψ when imposing derivatives at boundary
        Ψfunc(x, y) = @. sin(π*(x-1)/2) * sin(π*(y-1)/2) 
        DΨxfunc(x, y) = @. π/2 * cos(π*(x-1)/2) * sin(π*(y-1)/2) 
        DΨyfunc(x, y) = @. π/2 * sin(π*(x-1)/2) * cos(π*(y-1)/2) 

        n = (8, 8)
        reynolds = 0

        mesh = SpectralMesh2D(n)
        probl = Cavity4Sided(mesh, reynolds)

        Ψexact = [Ψfunc(x,y) for x in mesh.xnodes, y in mesh.ynodes] 

        Ψbcleft(y) = DΨxfunc(1, y)
        Ψbcright(y) = DΨxfunc(-1, y)
        Ψbctop(x) = DΨyfunc(x, 1)
        Ψbcbottom(x) = DΨyfunc(x, -1)

        setNeumannBC2D(probl, Ψbcleft, Ψbcright, Ψbctop, Ψbcbottom)

        Ψi = Ψexact[3:n[1]-1, 3:n[2]-1]
        Ψ = NS2DBenchmarkSolver.constructΨ(probl, Ψi)

        @test Ψ ≈ Ψexact atol=1e-6

        # Test right-hand-side function of equation for streamfunction 
        # in cavity flow
        n = (6, 6)
        mesh = SpectralMesh2D(n)

        reynolds = 100
        probl = Cavity4Sided(mesh, reynolds)

        k0 = 10
        bcfunc(x) = ((exp(k0*(x-1)) - 1) * (exp(-k0*(x+1)) - 1))^2

        bcleft = bcfunc.(probl.mesh.ynodes)
        bcright = -bcfunc.(probl.mesh.ynodes)
        bctop = bcfunc.(probl.mesh.xnodes)
        bcbottom = -bcfunc.(probl.mesh.xnodes)

        setNeumannBC2D(probl, bcleft, bcright, bctop, bcbottom)

        Ψi = zeros((n[1]-3, n[2]-3))

        Ψ = NS2DBenchmarkSolver.constructΨ(probl, Ψi)

        FΨ = NS2DBenchmarkSolver.rhs(probl, Ψi)
        FΨref = [-0.0624820681282460, 0.576027334953858, 0.245807608047586, -0.642691167651798, 0.0121206968541711, -0.642691167651798, 0.245807608047580, 0.576027334953857, -0.0624820681282459]
        @test FΨ ≈ FΨref

        # Test right-hand-side function with time stepping of equation for streamfunction 
        # in cavity flow
        ψi = ones((n[1]-3)*(n[2]-3))
        Ψold = ones((n[1]+1, n[2]+1))
        Δt = 1

        FΨ = NS2DBenchmarkSolver.rhstime(probl, Δt, Ψold, ψi)
        FΨref = [21.7652385080886, 7.16533446461262, 17.3663586367691, 3.45590293725733, -6.84022498215817, 3.45590293725730, 17.3663586367691, 7.16533446461265, 21.7652385080885]
        @test FΨ ≈ FΨref

        # Test to construct initial guess 
        Ψinitial = NS2DBenchmarkSolver.calculateΨinitial(probl, nbtimesteps=200)
        Ψinitialref = [0 0 0 0 0 0 0
                       0 0.0223534291102081 0.0655941971343996 0.0589510766581025 0.0555732229353761 0.0204248917264253 0
                       0 -0.0664329942047326 -0.0337131774973950 -0.0810657532138882 -0.0857836868595323 -0.0764539684037562 0
                       0 -0.0892424343792270 -0.150714381583066 -0.187874748844775 -0.150714381583067 -0.0892424343792272 0
                       0 -0.0764539684037564 -0.0857836868595334 -0.0810657532138900 -0.0337131774973966 -0.0664329942047329 0
                       0 0.0204248917264252 0.0555732229353758 0.0589510766581022 0.0655941971343994 0.0223534291102081 0
                       0 0 0 0 0 0 0]
        @test Ψinitial ≈ Ψinitialref
    end

    @testset "Complete setup and solve" begin
        n = (4, 4)
        nodesref = [1; 0.707106781186548; 0; -0.707106781186548; -1]
        D1ref = [5.5 -6.82842712474619 2  -1.17157287525381 0.5
                 1.70710678118655 -0.707106781186548 -1.41421356237310 0.707106781186548 -0.292893218813453
                 -0.5 1.41421356237310 0 -1.41421356237310 0.5
                 0.292893218813453 -0.707106781186548 1.41421356237310 0.707106781186548 -1.70710678118655
                 -0.5 1.17157287525381 -2 6.82842712474619 -5.5]

        mesh = SpectralMesh2D(n)

        reynolds = 500
        probl = Cavity4Sided(mesh, reynolds)

        @test probl.mesh.xnodes ≈ nodesref
        @test probl.mesh.ynodes ≈ nodesref
        @test probl.mesh.diffx1 ≈ D1ref
        @test probl.mesh.diffy1 ≈ D1ref

        # Test boundary construction 
        k0 = 10
        bcfunc(x) = @. ((exp(k0*(x-1)) - 1) * (exp(-k0*(x+1)) - 1))^2
        bcfuncneg(x) = .-bcfunc(x)

        setNeumannBC2D(probl, bcfuncneg, bcfunc, bcfunc, bcfuncneg)

        Ψref =  [0 0 0 0 0
                 0 0.0807493117423042 0.124977301580937 0.0807493117423042 0
                 0 -0.124977301580937 0 -0.124977301580937 0
                 0 0.0807493117423043 0.124977301580937 0.0807493117423042 0
                 0 0 0 0 0]

        Ψi = zeros((probl.mesh.nx-3, probl.mesh.ny-3))
        Ψ = NS2DBenchmarkSolver.constructΨ(probl, Ψi)
        @test Ψ ≈ Ψref

        # Test solve 
        sol = solve(probl)

        Ψsolref = [0 0 0 0 0
                   0 0.0726743805680738 0.0926775768840154 0.0726743805680738 0
                   0 -0.157277026277859 -0.129198898787687 -0.157277026277859 0
                   0 0.0726743805680738 0.0926775768840155 0.0726743805680738 0
                   0 0 0 0 0]

        @test sol.iter == 3
        @test sol.vals ≈ Ψsolref
    end
end
