function continuation_arclength(folder, Ψi, p::CavityStruct, Re_start, ΔRe, steps;
    save_steps = 1)
    @unpack n, scl = p.params

    # Create directories
    folderpsis = "$(folder)/psis"
    mkdir(folderpsis)

    # Write header
    fileresults = "$(folder)/results.csv"
    write_header(fileresults)

    p.params.Re = Re_start
    dim = (n - 3) * (n - 3)
    @inbounds u0 = reshape(Ψi[3:(n - 1), 3:(n - 1)], dim)
    u1, _, _ = newton(f!, u0, p)
    Ψ1 = constructBC(u1, p)
    saveΨ(folderpsis, Ψ1, 0, p.params.Re)
    save_result(fileresults, Ψ1, 0, 0, 0.0, p)

    p.params.Re = Re_start + ΔRe
    u2, _, _ = newton(f!, u0, p)
    Ψ2 = constructBC(u2, p)
    saveΨ(folderpsis, Ψ2, 1, p.params.Re)
    save_result(fileresults, Ψ2, 1, 0, 0.0, p)

    # Augmented system
    x1 = [u1; p.params.Re / scl]
    x2 = [u2; p.params.Re / scl]

    # Step size norm
    s = norm(x2 - x1)

    # Setup continuation cache
    x = similar(x1)
    fx = similar(x1)
    dx = similar(x1)
    v = similar(x1)
    xp = similar(x1)

    dim = (n - 3) * (n - 3)
    xi = zeros(dim)
    fxi = similar(xi)
    fxi_Re = similar(xi)

    J = zeros(dim, dim)
    J2 = zeros(dim + 1, dim + 1)
    jac_cache = FiniteDiff.JacobianCache(xi)

    cont_cache = (x, fx, dx, v, xp, xi, fxi, fxi_Re, J, J2, jac_cache)

    for i in 2:(steps)
        time = @elapsed tmp, iter, _ = newton_continuation(f!, x1, x2, s, p, cont_cache)

        x1 .= x2
        x2 .= tmp

        p.params.Re = tmp[end] * scl
        Ψ = constructBC(tmp[1:(end - 1)], p)

        if isinteger(i / save_steps)
            saveΨ(folderpsis, Ψ, i, p.params.Re)
        end
        save_result(fileresults, Ψ, i, iter, time, p)
    end
    return
end

function continuation_arclength_lsa(folderlsa, name, Ψi1, Ψi2, p::CavityStruct, Re1, Re2,
    steps; Re_stop_mode = 0)
    @unpack n, scl = p.params

    # Write header
    filelsa = "$folderlsa/results_$(name).csv"
    write_header_lsa(filelsa)

    p.params.Re = Re1
    @inbounds u1 = reshape(Ψi1[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    lsa_time = @elapsed lambdas = linearstability_lambdas(u1, p)
    save_result_lsa(filelsa, Ψi1, 0, lambdas, lsa_time, 0, 0, p)

    p.params.Re = Re2
    @inbounds u2 = reshape(Ψi2[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    lsa_time = @elapsed lambdas = linearstability_lambdas(u2, p)
    save_result_lsa(filelsa, Ψi2, 1, lambdas, lsa_time, 0, 0, p)

    # Augmented system
    x1 = [u1; p.params.Re / scl]
    x2 = [u2; p.params.Re / scl]

    # Step size norm
    s = norm(x2 - x1)

    # Setup continuation cache
    x = similar(x1)
    fx = similar(x1)
    dx = similar(x1)
    v = similar(x1)
    xp = similar(x1)

    dim = (n - 3) * (n - 3)
    xi = zeros(dim)
    fxi = similar(xi)
    fxi_Re = similar(xi)

    J = zeros(dim, dim)
    J2 = zeros(dim + 1, dim + 1)
    jac_cache = FiniteDiff.JacobianCache(xi)

    cont_cache = (x, fx, dx, v, xp, xi, fxi, fxi_Re, J, J2, jac_cache)

    for i in 2:(steps)
        newton_time = @elapsed tmp, iter, _ = newton_continuation(f!, x1, x2, s, p, cont_cache)

        x1 .= x2
        x2 .= tmp

        # Break after reaching again same Reynolds number
        # decreasing Re (=1) or increasing (=2)
        if Re_stop_mode == 1
            if tmp[end] * scl < Re1
                @warn("LSA stopped after $i of $steps steps: Reynolds number reached same value again")
                break
            end
        elseif Re_stop_mode == 2
            if tmp[end] * scl > Re1
                @warn("LSA stopped after $i of $steps steps: Reynolds number reached same value again")
                break
            end
        end

        p.params.Re = tmp[end] * scl
        Ψ = constructBC(tmp[1:(end - 1)], p)

        lsa_time = @elapsed lambdas = linearstability_lambdas(tmp[1:(end - 1)], p)
        save_result_lsa(filelsa, Ψ, i, lambdas, lsa_time, iter, newton_time, p)
    end
    return
end

function continuation_arclength_lsa(folderlsa, name, Ψi, p::CavityStruct, Re_start, ΔRe,
    steps)
    @unpack n = p.params

    Re1 = Re_start
    p.params.Re = Re1
    @inbounds u0 = reshape(Ψi[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    u1, _, _ = newton(f!, u0, p)
    Ψi1 = constructBC(u1, p)

    Re2 = Re_start + ΔRe
    p.params.Re = Re2
    u2, _, _ = newton(f!, u0, p)
    Ψi2 = constructBC(u2, p)

    return continuation_arclength_lsa(folderlsa, name, Ψi1, Ψi2, p, Re1, Re2, steps)
end

function newton_continuation(f!, x1, x2, s, p::CavityStruct, cont_cache;
    abstol = 1e-10, maxiters = 100)
    @unpack scl = p.params

    x, fx, dx, v, xp, xi, fxi, fxi_Re, J, J2, jac_cache = cont_cache

    @. x = x2
    @views @. xi = x[1:(end - 1)]

    dim = size(x, 1) - 1

    @. v = x2 - x1
    xp = x2 .+ (s / norm(v)) .* v

    # Fixed epsilon
    eps = 1e-8

    iter = 0
    tol = 1.0

    while tol > abstol && iter < maxiters
        p.params.Re = x[end] * scl # unscaled reynolds

        @views @. xi = x[1:(end - 1)]
        @views @. fxi = fx[1:(end - 1)]

        f!(fxi, xi, p)
        fx[1:(end - 1)] = fxi
        @. dx = x - xp
        fx[end] = dot(v, dx)

        FiniteDiff.finite_difference_jacobian!(J,
            (fxi, xi) -> f!(fxi, xi, p),
            xi,
            jac_cache)

        J2[1:dim, 1:dim] = J

        # Calculate change in "Reynolds variable"
        p.params.Re = (x[end] + eps) * scl # unscaled Reynolds
        f!(fxi_Re, xi, p)
        J2[1:dim, dim + 1] = (fxi_Re - fxi) / eps

        J2[dim + 1, :] = vec(v)

        dx = J2 \ fx
        @. x -= dx

        tol = norm(dx)
        iter += 1
    end

    if tol > abstol && iter == maxiters
        @warn("Newton for continuation did not converge in $iter iterations.")
    end

    return x, iter, tol
end

function newton_continuation(f!, x1, x2, s, p::CavityStruct; abstol = 1e-10, maxiters = 100)
    @unpack n = p.params

    x = similar(x1)
    fx = similar(x1)
    dx = similar(x1)
    v = similar(x1)
    xp = similar(x1)

    dim = (n - 3) * (n - 3)
    xi = zeros(dim)
    fxi = similar(xi)
    fxi_Re = similar(xi)

    J = zeros(dim, dim)
    J2 = zeros(dim + 1, dim + 1)
    jac_cache = FiniteDiff.JacobianCache(xi)

    cont_cache = (x, fx, dx, v, xp, xi, fxi, fxi_Re, J, J2, jac_cache)

    x, iter, tol = newton_continuation(f!,
        x1,
        x2,
        s,
        p,
        cont_cache;
        abstol = abstol,
        maxiters = maxiters)

    return x, iter, tol
end

function saveΨ(folder, Ψ, step, Re)
    writedlm("$folder/psi_step$(@sprintf("%03d", step))_Re$(@sprintf("%07.3f", Re)).txt",
        Ψ)

    return nothing
end

function save_result(fileresults, Ψ, step, iter, time, p)
    @unpack Re, D1, ic, i1, i2 = p.params

    U = D1 * Ψ
    V = Ψ * D1'

    # Important: indices are transposed, mapping to physical space
    result = "$step,$Re,$(Ψ[ic,ic]),$(Ψ[ic,i2]),$(Ψ[ic,i1]),$(Ψ[i2,ic])," *
             "$(U[ic,i2]),$(V[ic,i2]),$(U[ic,i1]),$(V[ic,i1]),$(U[i2,ic]),$(V[i2,ic])," *
             "$iter,$(time)\n"

    io = open(fileresults, "a")
    write(io, result)
    close(io)

    return nothing
end

function save_result_lsa(filelsa, Ψ, step, lambdas, lsa_time, iter, newton_time, p)
    @unpack Re, D1, ic, i1, i2 = p.params

    U = D1 * Ψ
    V = Ψ * D1'

    # Important: indices are transposed, mapping to physical space
    result = "$step,$Re,$(real(lambdas[1])),$(real(lambdas[2])),$(real(lambdas[3])),$(real(lambdas[4]))," *
             "$(real(lambdas[5])),$(real(lambdas[6])),$(real(lambdas[7])),$(real(lambdas[8]))," *
             "$(imag(lambdas[1])),$(imag(lambdas[2])),$(imag(lambdas[3])),$(imag(lambdas[4]))," *
             "$(imag(lambdas[5])),$(imag(lambdas[6])),$(imag(lambdas[7])),$(imag(lambdas[8]))," *
             "$(Ψ[ic,ic]),$(Ψ[ic,i2]),$(Ψ[ic,i1]),$(Ψ[i2,ic])," *
             "$(U[ic,i2]),$(V[ic,i2]),$(U[ic,i1]),$(V[ic,i1]),$(U[i2,ic]),$(V[i2,ic])," *
             "$lsa_time,$iter,$(newton_time)\n"

    io = open(filelsa, "a")
    write(io, result)
    close(io)

    return nothing
end

function write_header(fileresults)
    header = [
        "step",
        "Re",
        "psi_c",
        "psi_t",
        "psi_b",
        "psi_l",
        "u_t",
        "v_t",
        "u_b",
        "v_b",
        "u_l",
        "v_l",
        "newton_iters",
        "newton_time",
    ]
    writedlm(fileresults, reshape(header, 1, length(header)), ',')

    return nothing
end

function write_header_lsa(filelsa)
    header = [
        "step",
        "Re",
        "lambda1re",
        "lambda2re",
        "lambda3re",
        "lambda4re",
        "lambda5re",
        "lambda6re",
        "lambda7re",
        "lambda8re",
        "lambda1im",
        "lambda2im",
        "lambda3im",
        "lambda4im",
        "lambda5im",
        "lambda6im",
        "lambda7im",
        "lambda8im",
        "psi_c",
        "psi_t",
        "psi_b",
        "psi_l",
        "u_t",
        "v_t",
        "u_b",
        "v_b",
        "u_l",
        "v_l",
        "lsa_time",
        "newton_iters",
        "newton_time",
    ]
    writedlm(filelsa, reshape(header, 1, length(header)), ',')

    return nothing
end
