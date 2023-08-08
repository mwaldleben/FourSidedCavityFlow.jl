"""
    f_linearstability!(fu, Ψ, Ψ0, p)

Computes the for a given solution Ψ0 and the eigenfunction Ψ. This function calculates
the elements of the matrix ``A`` generalized eigenvalue problem. Stores the solution  
in `fu`.
"""
function f_linearstability!(fu, Ψ, Ψ0, p::CavityStruct)
    @unpack Re, n, D1, D2, D4 = p.params
    @unpack fΨ, D2Ψ, ΨD2, D4Ψ, ΨD4, laplΨ, biharmΨ = p.cache

    # Ψ0
    D1Ψ0 = similar(Ψ)
    Ψ0D1 = similar(Ψ)
    D2Ψ0 = similar(Ψ)
    Ψ0D2 = similar(Ψ)
    laplΨ0 = similar(Ψ)
    D1laplΨ0 = similar(Ψ)
    laplΨ0D1 = similar(Ψ)

    mul!(D1Ψ0, D1, Ψ0)
    mul!(Ψ0D1, Ψ0, D1')

    mul!(D2Ψ0, D2, Ψ0)
    mul!(Ψ0D2, Ψ0, D2')
    @inbounds @. laplΨ0 = D2Ψ0 + Ψ0D2

    mul!(D1laplΨ0, D1, laplΨ0)
    mul!(laplΨ0D1, laplΨ0, D1')

    # Ψ
    D1Ψ = similar(Ψ)
    ΨD1 = similar(Ψ)
    D1laplΨ = similar(Ψ)
    laplΨD1 = similar(Ψ)

    mul!(D1Ψ, D1, Ψ)
    mul!(ΨD1, Ψ, D1')

    mul!(D2Ψ, D2, Ψ)
    mul!(ΨD2, Ψ, D2')
    @inbounds @. laplΨ = D2Ψ + ΨD2

    mul!(D1laplΨ, D1, laplΨ)
    mul!(laplΨD1, laplΨ, D1')

    mul!(D4Ψ, D4, Ψ)
    mul!(ΨD4, Ψ, D4')
    mul!(laplΨ, D2Ψ, D2') # using as intermediate memory
    @inbounds @. biharmΨ = D4Ψ + ΨD4 + 2 * laplΨ

    @inbounds @. fΨ = (1 / Re) * biharmΨ + D1Ψ0 .* laplΨD1 - Ψ0D1 .* D1laplΨ +
                      D1Ψ .* laplΨ0D1 - ΨD1 .* D1laplΨ0

    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]

    return nothing
end

"""
    f_linearstability!(fu, Ψ, Ψ0, p)

Computes the matrices `A` and `B` of the linear stability around a base state `u`.  
"""
function linearstability_matrices!(A, B, u, p::CavityStruct)
    @unpack Re, n, D2 = p.params
    @unpack D2Ψ, ΨD2, laplΨ = p.cache

    Ψ0 = constructBC(u, p)

    fu = similar(u)
    for i in 1:(n - 3)
        for j in 1:(n - 3)
            Ψ = zeros(n + 1, n + 1)
            Ψ[i + 2, j + 2] = 1.0

            construct_homogenousBC!(Ψ, p)
            f_linearstability!(fu, Ψ, Ψ0, p)
            k = (j - 1) * (n - 3) + i
            @views A[:, k] .= fu

            mul!(D2Ψ, D2, Ψ)
            mul!(ΨD2, Ψ, D2')
            @inbounds @. laplΨ = D2Ψ + ΨD2
            @views B[:, k] .= laplΨ[3:(n - 1), 3:(n - 1)][:]
        end
    end

    return nothing
end

"""
    linearstability_lambdas(u, p)

Computes the eigenvalues of the linear stability around a base state `u`. 
"""
function linearstability_lambdas(u, p::CavityStruct)
    @unpack Re, n = p.params

    dim = (n - 3) * (n - 3)
    A = zeros(dim, dim)
    B = zeros(dim, dim)

    linearstability_matrices!(A, B, u, p)

    vals, _ = eigen(A, B) # generalized eigenvalues

    # Sorting complex lambdas lexicographically
    lambdas = sort(vals; by = λ -> (floor(real(λ); digits = 5), imag(λ)), rev = true)
    return lambdas
end

"""
    linearstability_lambdamax(Re, u, p)

Computes the real part of the rightmost eigenvalues of the linear stability around a base state `u` and a Reynolds number `Re`. 
"""
function linearstability_lambdamax(Re, u, p::CavityStruct)
    @unpack n = p.params

    dim = (n - 3) * (n - 3)
    A = zeros(dim, dim)
    B = zeros(dim, dim)

    p.params.Re = Re
    linearstability_matrices!(A, B, u, p)

    vals, _ = eigen(A, B) # generalized eigenvalues
    λmax = sort(real(vals); rev = true)[1]

    return λmax
end

"""
    Re, u, iter, tol = newton1D_for_linearstability(Re0, u0, p; abstol = 1e-10, maxiters = 20, verbose = false)

    1D Newton to calculates the Reynolds number `Re` where the rightmost eigenvalues real part is 0. `u0` is the initial guess where
    for the bifurcation. It only works if the state where the bifurcation is has a well-conditioned Jacobian.    
"""
function newton1D_for_linearstability(Re0, u0, p::CavityStruct; abstol = 1e-10,
    maxiters = 20, verbose = false)
    @unpack n, scl = p.params

    Re = Re0
    x = Re0 / scl

    u = copy(u0)

    iter = 0
    tol = 1.0

    if verbose == true
        @printf("  %-10s %-10s %-10s %-10s\n", "Newtonstep", "Re", "lambda_max", "Time[s]")
    end

    while tol > abstol && iter < maxiters
        time = @elapsed begin
            # Refine solution for new Reynolds number
            p.params.Re = Re
            u, _, _ = newton(f!, u, p)
            fx = linearstability_lambdamax(Re, u, p)

            # Fixed step size!
            x1 = x + 1e-8
            Re1 = x1 * scl

            # Refine solution for step
            p.params.Re = Re1
            u1, _, _ = newton(f!, u, p)
            fx1 = linearstability_lambdamax(Re1, u1, p)

            dfx = (fx1 - fx) / 1e-8

            dx = -fx / dfx
            x = x + dx
            Re = x * scl

            tol = abs(dx)
            iter += 1
        end

        p.params.Re = Re
        u, _, _ = newton(f!, u, p; abstol = 1e-9, maxiters = 150)

        if verbose == true
            @printf("  %-10d %-10.6f %-10.6f %-10.6f\n", iter, Re, fx1, time)
        end
    end

    return Re, u, iter, tol
end
