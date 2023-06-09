function continuation_arclength(folder, Ψi, p::CavityStruct, Re_start, ΔRe, steps; save_steps = 1)
    @unpack n, scl = p.params

    # Create directories
    folderpsis = "$(folder)/psis"
    mkdir(folderpsis)

    # Write header and open results file
    fileresults = "$(folder)/results.csv"
    header = ["step", "Re","psi_c", "psi_t", "psi_b", "psi_l",
              "u_t", "v_t", "u_b", "v_b", "u_l", "v_l",
              "newton_iters", "newton_time"]
    writedlm(fileresults, reshape(header,1,length(header)), ',')
    io = open(fileresults, "a")

    p.params.Re = Re_start
    @inbounds u0 = reshape(Ψi[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    u1, _, _ = newton(f!, u0, p)
    Ψ1 = constructBC(u1, p)
    saveΨ(folderpsis, Ψ1, 0, p.params.Re)
    save_result(io, Ψ1, 0, 0, 0.0, p)
    u1 = [u1; p.params.Re / scl]

    p.params.Re = Re_start + ΔRe
    u2, _, _ = newton(f!, u0, p)
    Ψ2 = constructBC(u2, p)
    saveΨ(folderpsis, Ψ2, 1, p.params.Re)
    save_result(io, Ψ2, 1, 0, 0.0, p)
    u2 = [u2; p.params.Re / scl]

    # Step size norm
    s = norm(u2 - u1)

    for i in 2:(steps)
        time = @elapsed u, iter, _ = newton_for_continuation(f!, u1, u2, s, p)

        u1 = u2
        u2 = u

        p.params.Re = u[end] * scl
        Ψ = constructBC(u[1:end-1], p)

        if isinteger(i / save_steps)
            saveΨ(folderpsis, Ψ, i, p.params.Re)
        end
        save_result(io, Ψ, i, iter, time, p)
    end
    close(io)
end

function saveΨ(folder, Ψ, step, Re)
    writedlm("$folder/psi_step$(@sprintf("%03d", step))_Re$(@sprintf("%07.3f", Re)).txt",  Ψ)
end

function save_result(io, Ψ, step, iter, time, p)
    @unpack Re, D1, ic, i1, i2 = p.params 

    U = D1*Ψ 
    V = Ψ*D1' 

    # Important: indices are transposed, mapping to physical space
    result = "$step,$Re,$(Ψ[ic,ic]),$(Ψ[ic,i2]),$(Ψ[ic,i1]),$(Ψ[i2,ic])," *
             "$(U[ic,i2]),$(V[ic,i2]),$(U[ic,i1]),$(V[ic,i1]),$(U[i2,ic]),$(V[i2,ic])," *
             "$iter,$(time)\n"
    write(io, result)
    flush(io)
end

function continuation_arclength_lsa(folderlsa, name, Ψi1, Ψi2, p::CavityStruct, Re1, Re2, steps;
                                    Re_stop_mode = 0)
    @unpack n, scl = p.params

    # Write header and open results file
    filelsa = "$folderlsa/results_$(name).csv"
    header = ["step", "Re", "lambda1", "lambda2", "lambda3", "lambda4", "lambda5", "lambda6", "lambda7", "lambda8",
              "psi_c", "psi_t", "psi_b", "psi_l",
              "u_t", "v_t", "u_b", "v_b", "u_l", "v_l",
              "lsa_time", "newton_iters", "newton_time"]
    writedlm(filelsa, reshape(header,1,length(header)), ',')
    io = open(filelsa, "a")

    p.params.Re = Re1
    @inbounds u1 = reshape(Ψi1[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    u1 = [u1; Re1 / scl]
    lsa_time = @elapsed lambdas = linearstability_lambdas(u1[1:end-1], p)
    save_result_lsa(io, Ψi1, 0, lambdas, lsa_time, 0, 0, p)

    p.params.Re = Re2
    @inbounds u2 = reshape(Ψi2[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    u2 = [u2; Re2 / scl]
    lsa_time = @elapsed lambdas = linearstability_lambdas(u2[1:end-1], p)
    save_result_lsa(io, Ψi2, 1, lambdas, lsa_time, 0, 0, p)

    # Step size norm
    s = norm(u2 - u1)

    for i in 2:(steps)
        newton_time = @elapsed u, iter, _ = newton_for_continuation(f!, u1, u2, s, p)

        u1 = u2
        u2 = u

        # Break after decreasing Reynoldsnumber (=1)
        # or increasing (=2)
        if Re_stop_mode == 1
            if u2[end] < u1[end] 
                @warn("LSA stopped after $i of $steps steps: Reynolds number is decreasing")
                break
            end
        elseif Re_stop_mode == 2
            if u2[end] > u1[end] 
                @warn("LSA stopped after $i of $steps steps: Reynolds number is increasing")
                break
            end
        end

        p.params.Re = u[end] * scl
        Ψ = constructBC(u[1:end-1], p)

        lsa_time = @elapsed lambdas = linearstability_lambdas(u[1:end-1], p)
        save_result_lsa(io, Ψ, i, lambdas, lsa_time, iter, newton_time, p)
    end
    close(io)
end

function continuation_arclength_lsa(folderlsa, name, Ψi, p::CavityStruct, Re_start, ΔRe, steps)
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

    continuation_arclength_lsa(folderlsa, name, Ψi1, Ψi2, p, Re1, Re2, steps)
end

function save_result_lsa(io, Ψ, step, lambdas, lsa_time, iter, newton_time, p)
    @unpack Re, D1, ic, i1, i2 = p.params 

    U = D1*Ψ 
    V = Ψ*D1' 

    # Important: indices are transposed, mapping to physical space
    result = "$step,$Re,$(lambdas[1]),$(lambdas[2]),$(lambdas[3]),$(lambdas[4])," *
             "$(lambdas[5]),$(lambdas[6]),$(lambdas[8]),$(lambdas[9])," *
             "$(Ψ[ic,ic]),$(Ψ[ic,i2]),$(Ψ[ic,i1]),$(Ψ[i2,ic])," *
             "$(U[ic,i2]),$(V[ic,i2]),$(U[ic,i1]),$(V[ic,i1]),$(U[i2,ic]),$(V[i2,ic])," *
             "$lsa_time,$iter,$(newton_time)\n"
    write(io, result)
    flush(io)
end

function newton_for_continuation(f!, x1, x2, s, p::CavityStruct; tolmax = 1e-10,
                                 maxiter = 100)
    @unpack scl = p.params
    x = copy(x2)
    xi = x[1:(end - 1)]

    dim = size(x, 1) - 1

    eps = 1e-8

    v = x2 - x1
    xp = x2 + (s / norm(v)) * v

    fx = similar(x)
    fxi_Re = zeros(dim)

    J = zeros(dim, dim)
    J2 = zeros(dim + 1, dim + 1)

    dx = zeros(dim + 1)
    cache = FiniteDiff.JacobianCache(xi)

    iter = 0
    tol = 1.0
    while tol > tolmax && iter < maxiter
        p.params.Re = x[end] * scl # unscaled reynolds

        xi = x[1:(end - 1)]
        fxi = fx[1:(end - 1)]

        f!(fxi, xi, p)
        fx[1:(end - 1)] = fxi
        fx[end] = v' * (x - xp)

        FiniteDiff.finite_difference_jacobian!(J, (fxi, xi) -> f!(fxi, xi, p), xi, cache)
        # jacobian!(J, f!, xi, p; dx = eps)

        J2[1:dim, 1:dim] = J

        # Calculate change in "Reynolds variable"
        p.params.Re = (x[end] + eps) * scl # unscaled Reynolds
        f!(fxi_Re, xi, p)
        J2[1:dim, dim + 1] = (fxi_Re - fxi) / eps

        J2[dim + 1, :] = vec(v)

        dx = J2 \ (-fx)
        @. x = x + dx

        tol = norm(dx)
        iter += 1
    end

    return x, iter, tol
end
