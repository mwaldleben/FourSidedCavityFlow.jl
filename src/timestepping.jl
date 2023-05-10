function timestepping(Ψstart, p::CavityParameters, Δt, nb_timesteps)
    @unpack n, Ψ, Ψ0 = p

    @inbounds Ψ0 .= Ψstart
    @inbounds u0 = reshape(Ψstart[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    ft!(fu, u, p) = ftime!(fu, u, p, Δt)

    for step in 1:nb_timesteps
        u0, iter, tol = newton(ft!, u0, p)

        Ψ[3:(n - 1), 3:(n - 1)] .= reshape(u0, (n - 3, n - 3))

        Ψ0 .= Ψ
    end

    construct_BC!(p)

    return Ψ
end

function timestepping_save(Ψstart, p::CavityParameters, Δt, steps)
    @unpack n, Ψ, Ψ0 = p

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
        construct_BC!(p)

        sol[i + 1] = Ψ
        time_series[i + 1] = Δt * i
    end

    return sol, time_series
end
