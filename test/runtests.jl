using NS2DBenchmarkSolver
using Test

@testset "differentation.jl" begin
    n = 2
    Dₓ = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
    x = [-1; 0; 1]
    @test NS2DBenchmarkSolver.diffchebychev(n)[1] ≈ x
    @test NS2DBenchmarkSolver.diffchebychev(n)[2] ≈ Dₓ
end

@testset "boundarycondition.jl" begin
    val = 1.1
    bc1 = BCDirichlet1D(val)
    bc2 = BCNeumann1D(val)
    @test bc1.val == val
    @test bc2.val == val
end

@testset "solution.jl" begin
    @testset "Solution1D" begin 
        x = [0; 1; 2]
        y = [1; 1; 1]

        sol = Solution1D(x, y)

        @test sol.nodes == x
        @test sol.vals == y 
    end
    @testset "Solution1D plot recipes" begin
        @test skip=true
    end
end

@testset "spectralmesh.jl" begin
    @testset "SpectralMesh1D" begin
        nbcells = 2
        span = (-1, 1)
        Dₓ = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
        x = [-1; 0; 1]

        mesh = SpectralMesh1D(nbcells, span)

        @test mesh.nodes ≈ x
        @test mesh.diff1mat ≈ Dₓ
        @test mesh.spectralmethod == "chebychev"
    end
    @testset "SpectralMesh2D" begin
        nbcells = (2, 2)
        xspan = (-1, 1)
        yspan = (-1, 1)
        Dₓ = [-1.5 2.0 -0.5; -0.5 0.0 0.5; 0.5 -2.0 1.5]
        x = [-1; 0; 1]

        mesh = SpectralMesh2D(nbcells, xspan, yspan)

        @test mesh.xnodes ≈ x
        @test mesh.ynodes ≈ x
        @test mesh.diffx1mat ≈ Dₓ
        @test mesh.diffy1mat ≈ Dₓ
        @test mesh.spectralmethod == "chebychev"
    end
end
@testset "boundarycondition.jl" begin
end

@testset "spectralproblem.jl" begin
    @testset "ExamplePoisson1D" begin
        nbcells = 16
        span = (-1, 1)

        mesh = SpectralMesh1D(nbcells, span)

        bcmin = BCDirichlet1D(0.0) 
        bcmax = BCDirichlet1D(0.0) 

        func(x) = exp(4*x) 
        probl = NS2DBenchmarkSolver.ExamplePoisson1D(mesh, func, (bcmin, bcmax)) 

        @test probl.bcmin.val == bcmin.val
        @test probl.bcmax.val == bcmax.val

        sol = NS2DBenchmarkSolver.solve(probl)

        exactsol(x) = @. (exp(4*x) - sinh(4)*x - cosh(4))/16
        exactvals = exactsol(mesh.nodes)

        @test sol.vals ≈ exactvals
    end
    @testset "Cavity4Sided functions" begin
        # Test constructBChelpermat function
        #nodes = [1; 1/2; -1/2; -1] 
        D = [19/6 -4 4/3 -1/2; 1 -1/3 -1 1/3; -1/3 1 1/3 -1; 1/2 -4/3 4 -19/6]

        Minvref = [-9 3; -3 9]/32
        Minv = NS2DBenchmarkSolver.constructBChelpermat(D)

        @test Minv ≈ Minvref

        # Test fill interior of Ψ function
        Ψ = NS2DBenchmarkSolver.fillΨint(8, 8, initialΨval=1)

        @test Ψ[3:8-1,3:8-1] == ones(5, 5) 

        # Test boundary reconstruction Ψ when imposing derivatives at boundary
        Ψfunc(x, y) = @. sin(π*(x-1)/2) * sin(π*(y-1)/2) 
        DΨxfunc(x, y) = @. π/2 * cos(π*(x-1)/2) * sin(π*(y-1)/2) 
        DΨyfunc(x, y) = @. π/2 * sin(π*(x-1)/2) * cos(π*(y-1)/2) 

        n = 4
        mesh = SpectralMesh2D((n, n))
        Ψexact = [Ψfunc(x,y) for x in mesh.xnodes, y in mesh.ynodes] 

        Ψbcxmin(y) = DΨxfunc(-1, y)
        Ψbcxmax(y) = DΨxfunc(1, y)
        Ψbcymin(x) = DΨyfunc(x, -1)
        Ψbcymax(x) = DΨyfunc(x, 1)

        bcxmin = BCNeumann2D(Ψbcxmin(mesh.ynodes))
        bcxmax = BCNeumann2D(Ψbcxmax(mesh.ynodes))
        bcymin = BCNeumann2D(Ψbcymin(mesh.xnodes))
        bcymax = BCNeumann2D(Ψbcymax(mesh.xnodes))

        Ψ = zeros(n+1, n+1)  
        Ψ[3:n-1,3:n-1] = Ψexact[3:n-1,3:n-1] 

        Ψ = NS2DBenchmarkSolver.constructΨboundary(Ψ, mesh.diffx1mat, mesh.diffy1mat, bcxmin, bcxmax, bcymin, bcymax)

        @test Ψ ≈ Ψexact atol=0.01
    end

    @testset "Cavity4Sided" begin
        nbcells = (4, 4)
        xspan = (1, -1)
        yspan = (1, -1)

        nodes = [1; 0.707106781186548; 0; -0.707106781186548; -1]
        D = [5.5 -6.82842712474619 2	-1.17157287525381	0.5
             1.70710678118655	-0.707106781186548 -1.41421356237310 0.707106781186548 -0.292893218813453
             -0.5 1.41421356237310 0 -1.41421356237310 0.5
             0.292893218813453 -0.707106781186548	1.41421356237310 0.707106781186548 -1.70710678118655
             -0.5 1.17157287525381 -2	6.82842712474619 -5.5]

        mesh = SpectralMesh2D(nbcells, xspan, yspan)

        probl = Cavity4Sided(mesh)

        @test probl.mesh.xnodes ≈ nodes
        @test probl.mesh.ynodes ≈ nodes
        @test probl.mesh.diffx1mat ≈ D
        @test probl.mesh.diffy1mat ≈ D

        # Test boundary construction 
        k0 =10
        bcfunc(x) = ((exp(k0*(x-1)) - 1) * (exp(-k0*(x+1)) -1))^2

        bcxmin = BCNeumann2D(-bcfunc.(probl.mesh.ynodes))
        bcxmax = BCNeumann2D(bcfunc.(probl.mesh.ynodes))
        bcymin = BCNeumann2D(bcfunc.(probl.mesh.xnodes))
        bcymax = BCNeumann2D(-bcfunc.(probl.mesh.xnodes))

        applyBC2D(probl, bcxmin, bcxmax, bcymin, bcymax)

        Ψref =  [0 0 0 0 0
                 0 0.0807493117423042	0.124977301580937	0.0807493117423042 0
                 0 -0.124977301580937	0	-0.124977301580937 0
                 0 0.0807493117423043	0.124977301580937	0.0807493117423042 0
                 0 0 0 0 0]

        @test probl.psiinit ≈ Ψref
    end
end
