function timestepping(Ψinit, p::CavityStruct, h, timesteps; save = false, verbose = false)
    @unpack n, ic, Ψi = p.params
    @unpack Ψ, Ψ0 = p.cache

    @inbounds Ψi .= Ψinit
    @inbounds u = reshape(Ψinit[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    ft!(fu, u, p) = ftime!(fu, u, p, h)

    if verbose == true
        @printf("  %-10s %-10s %-10s %-10s\n", "Timestep", "Ψc", "Newton[s]", "Iters.")
    end

    if save == true
        sol = Vector{typeof(Ψi)}(undef, timesteps)
        time = Vector{Float64}(undef, timesteps)

        sol[1] = Ψi
        time[1] = 0
    end

    for ts in 1:timesteps
        newton_time = @elapsed u, iter, tol = newton(ft!, u, p)

        Ψ[3:(n - 1), 3:(n - 1)] .= reshape(u, (n - 3, n - 3))

        constructBC!(Ψ, p)

        if verbose == true
            @printf("  %-10d %-+10.6f %-10.6f %-10d\n", ts, Ψi[ic, ic], newton_time, iter)
        end

        if save == true
            sol[ts + 1] = copy(Ψ)
            time[ts + 1] = h * ts
        end

        Ψi .= Ψ
    end

    if save == true
        return sol, time
    else
        return Ψi
    end
end
