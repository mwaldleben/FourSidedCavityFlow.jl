"""
    Ψi = timestepping(Ψinit, p, h, timesteps; savesteps = false, verbose = false)

Timestepping of the equation of motion starting with an initial solution `Ψinit` and then perfoms a specified
number of timesteps using the implicit Euler scheme.

# Arguments
- `Ψinit::Matrix`: Initial streamfunction field to start time-integration 
- `p::CavityStruct`: Parameters and cache for the problem. 
- `h::Real`: Timestep size 
- `timesteps:Integer`: Number of timesteps
- `savesteps::Bool=false`: If false only returns the solution of the last timesteps, if true a solution vector of all timesteps and the timeseries
- `verbose::Bool=false`: If true prints information of the state at each timesteps
"""
function timestepping(Ψinit, p::CavityStruct, h, timesteps; savesteps = false, verbose = false)
    @unpack n, ic, Ψi = p.params
    @unpack Ψ = p.cache

    @inbounds Ψi .= Ψinit
    @inbounds u = reshape(Ψinit[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    ft!(fu, u, p) = ftime!(fu, u, p, h)

    if verbose == true
        @printf("  %-10s %-10s %-10s %-10s %-10s\n",
            "Timestep", "Time", "Ψc", "Newton[s]", "Iters.")
        @printf("  %-10d %-10.6f %-+10.6f %-10.6f %-10d\n", 0, 0.0, Ψi[ic, ic], 0.0, 0)
    end

    if savesteps == true
        sol = Array{Float64}(undef, (timesteps + 1, n + 1, n + 1))
        timeseries = Vector{Float64}(undef, timesteps + 1)

        sol[1, :, :] = Ψi
        timeseries[1] = 0
    end

    # Newton cache
    x = similar(u)
    fx = similar(u)
    dx = similar(u)

    dim = size(u, 1)
    J = zeros(dim, dim)
    jac_cache = FiniteDiff.JacobianCache(u)

    newton_cache = (x, fx, dx, J, jac_cache)

    for ts in 1:timesteps
        newton_time = @elapsed u, iter, tol = newton(ft!, u, p, newton_cache)

        Ψ[3:(n - 1), 3:(n - 1)] .= reshape(u, (n - 3, n - 3))

        constructBC!(Ψ, p)

        if savesteps == true
            sol[ts + 1, :, :] = copy(Ψ)
            timeseries[ts + 1] = h * ts 
        end

        if verbose == true
            @printf("  %-10d %-10.6f %-+10.6f %-10.6f %-10d\n",
                ts, h * ts, Ψ[ic, ic], newton_time, iter)
        end

        Ψi .= Ψ
    end

    if savesteps == true
        return sol, timeseries
    else
        return Ψi
    end
end
