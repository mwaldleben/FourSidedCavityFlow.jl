function timestepping(Ψstart, p::CavityStruct, Δt, timesteps; convergence_check=false, verbose=false)
    @unpack n = p.params
    @unpack Ψ, Ψ0 = p.cache

    @inbounds Ψ0 .= Ψstart
    @inbounds u0 = reshape(Ψstart[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    ft!(fu, u, p) = ftime!(fu, u, p, Δt)

    nc = Int(ceil((n-3)*(n-3)/2)) 

    if verbose == true
        @printf("  %-10s %-10s %-10s %-10s\n", "Timestep", "Ψc", "Newton[s]", "Iters.")
    end

    for timestep in 1:timesteps
        time = @elapsed u, iter, tol = newton(ft!, u0, p)

        Ψ[3:(n - 1), 3:(n - 1)] .= reshape(u, (n - 3, n - 3))

        if verbose == true
            @printf("  %-10d %-+10.6f %-10.6f %-10d\n", timestep, u[nc], time, iter)
        end

        if convergence_check == true
            if abs(u0[nc] - u[nc]) < 1e-8
                println(abs(u0[nc] - u[nc]))
                @printf("  ... center value converged")
                break
            end
        end

        u0 .= u

        # TODO: why error when removing
        Ψ0 .= Ψ 
    end

    constructBC!(Ψ, p)

    return Ψ
end

function timestepping_save(Ψstart, p::CavityStruct, Δt, steps)
    @unpack n = p.params
    @unpack n, Ψ, Ψ0 = p.cache

    @inbounds Ψ0 .= Ψstart
    @inbounds u0 = reshape(Ψstart[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))
    fu = similar(u0)

    ft!(fu, u, p) = ftime!(fu, u, p, Δt)

    sol = Vector{typeof(Ψstart)}(undef, steps)
    time_series = Vector{Float64}(undef, steps)

    sol[1] = Ψstart
    time_series[1] = 0

    for i in 1:(steps - 1)
        u0, iter, tol = newton(ft!, u0, p)

        Ψ[3:(n - 1), 3:(n - 1)] .= reshape(u0, (n - 3, n - 3))
        constructBC!(Ψ, p)

        sol[i + 1] = Ψ
        time_series[i + 1] = Δt * i
    end

    return sol, time_series
end
