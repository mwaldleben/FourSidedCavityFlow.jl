function f_linearstability!(fu, Ψ, Ψ0, p::CavityParameters)
    @unpack Re, n, D1, D2, D4, fΨ, D2Ψ, ΨD2, D4Ψ, ΨD4, laplΨ, biharmΨ, laplΨ = p

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

    # TODO: fix stream function formulation, problem with linebreak as well
    # @inbounds @. fΨ = (1 / Re) * biharmΨ + D1Ψ0 * laplΨD1 
    # + D1Ψ * laplΨ0D1 - Ψ0D1 * D1laplΨ - ΨD1 * D1laplΨ0 
    # @inbounds @. fΨ = (1 / Re) * biharmΨ - D1Ψ0 * laplΨD1 - D1Ψ * laplΨ0D1 + Ψ0D1 * D1laplΨ + ΨD1 * D1laplΨ0 
    @inbounds @. fΨ = (1 / Re) * biharmΨ - (D1Ψ0 .* laplΨD1 - Ψ0D1 .* D1laplΨ) -
                      (D1Ψ .* laplΨ0D1 - ΨD1 .* D1laplΨ0)

    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]

    return nothing
end

function linearstability_matrices!(A, B, u, p)
    @unpack Re, n, D2, D2Ψ, ΨD2, laplΨ = p

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

function linearstability_lambdas!(lambdas_real, u, p)
    @unpack Re, n = p

    dim = (n - 3) * (n - 3)
    A = zeros(dim, dim)
    B = zeros(dim, dim)

    linearstability_matrices!(A, B, u, p)

    vals, vecs = eigen(A, B) # generalized eigenvalues
    lambdas_real .= sort(real(vals), rev = true)
end

function linearstability_lambdas(u, p)
    @unpack Re, n = p

    dim = (n - 3) * (n - 3)
    A = zeros(dim, dim)
    B = zeros(dim, dim)

    linearstability_matrices!(A, B, u, p)

    vals, vecs = eigen(A, B) # generalized eigenvalues
    lambdas_real = sort(real(vals), rev = true)
    return lambdas_real
end

function linearstability_lambdamax(Re, u, p)
    @unpack n, Ψ = p

    dim = (n - 3) * (n - 3)
    A = zeros(dim, dim)
    B = zeros(dim, dim)

    p.Re = Re
    linearstability_matrices!(A, B, u, p)

    vals, _ = eigen(A, B) # generalized eigenvalues
    λmax = sort(real(vals), rev = true)[1]

    return λmax
end

function newton1D_for_linearstability(Re0, u0, p; tolmax = 1e-10, maxiter = 20, verbose = false)
    @unpack n, Ψ, scl = p

    Re = Re0
    x = Re0 / scl

    u = copy(u0)

    iter = 0
    tol = 1.0

    if verbose == true
        @printf("  %-10s %-10s %-10s\n", "Newtonstep", "Re", "Time[s]")
    end

    while tol > tolmax && iter < maxiter
        time = @elapsed begin 
            # Refine solution for new Reynolds number
            p.Re = Re
            u, _, _ = newton(f!, u, p)
            fx = linearstability_lambdamax(Re, u, p)

            # Fixed step size!
            x1 = x + 1e-8
            Re1 = x1 * scl

            # Refine solution for step
            p.Re = Re1
            u1, _, _ = newton(f!, u, p)
            fx1 = linearstability_lambdamax(Re1, u1, p)

            dfx = (fx1 - fx) / 1e-8

            dx = -fx / dfx
            x = x + dx
            Re = x * scl

            tol = abs(dx)
            iter += 1
        end

        if verbose == true
            @printf("  %-10d %-10.6f %-10.6f\n", iter, Re, time)
        end
    end

    return Re, iter, tol
end
