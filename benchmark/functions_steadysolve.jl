using LinearAlgebra
using UnPack
using FiniteDiff
using BenchmarkTools
using Random
using Plots;
gr();

mutable struct CavityCache{T}
    # Parameters
    fΨ::Matrix{T}
    Ψ::Matrix{T}
    n::Int64
    Re::T
    D1::Matrix{T}
    D2::Matrix{T}
    D4::Matrix{T}
    bcleft::Vector{T}
    bcright::Vector{T}
    bctop::Vector{T}
    bcbottom::Vector{T}

    # Cache
    Minv::Matrix{T}
    h1::Vector{T}
    h2::Vector{T}
    k1::Vector{T}
    k2::Vector{T}
    DΨ::Matrix{T}
    ΨD::Matrix{T}
    lDΨ::Matrix{T}
    lΨD::Matrix{T}
    biharmΨ::Matrix{T}
    laplΨ::Matrix{T}
    nonlinΨ::Matrix{T}
end

function diffchebyshev(n::Int)
    # Define nodes in [1, -1]
    nodes = [cos(k * π / n) for k in 0:n]

    # Off-diagonal entries
    c = [2; ones(n - 1); 2]
    dij = (i, j) -> (-1)^(i + j) * c[i + 1] / (c[j + 1] * (nodes[i + 1] - nodes[j + 1]))
    D1 = [dij(i, j) for i in 0:n, j in 0:n]

    # Diagonal entries
    D1[isinf.(D1)] .= 0
    s = sum(D1, dims = 2)
    D1 -= diagm(s[:, 1])

    # Second and fourth derivative
    D2 = D1^2
    D4 = D2^2

    return nodes, D1, D2, D4
end

function constructBCmatrix(D)
    n = size(D, 1) - 1
    M = [D[1, 2] D[1, n]; D[n + 1, 2] D[n + 1, n]]
    return inv(M)
end

function constructΨ!(Ψ, n, D, Minv, bcleft, bcright, bctop, bcbottom, h1, h2, k1, k2)
    @views Ψ1 = Ψ[3:(n - 1), 3:(n - 1)]
    @views Ψ2 = Ψ[3:(n - 1), 2:n]'
    @views Dcs = D[1, 3:(n - 1)]
    @views Dce = D[n + 1, 3:(n - 1)]

    mul!(h1, Ψ1, Dcs)
    mul!(h2, Ψ1, Dce)

    @views @. Ψ[3:(n - 1), 2] = Minv[1, 1] * (bctop[3:(n - 1)] - h1) +
                                Minv[1, 2] * (bcbottom[3:(n - 1)] - h2)
    @views @. Ψ[3:(n - 1), n] = Minv[2, 1] * (bctop[3:(n - 1)] - h1) +
                                Minv[2, 2] * (bcbottom[3:(n - 1)] - h2)

    mul!(k1, Ψ2, Dcs)
    mul!(k2, Ψ2, Dce)

    @views @. Ψ[2, 2:n] = Minv[1, 1] * (bcright[2:n] - k1) + Minv[1, 2] * (bcleft[2:n] - k2)
    @views @. Ψ[n, 2:n] = Minv[2, 1] * (bcright[2:n] - k1) + Minv[2, 2] * (bcleft[2:n] - k2)

    return nothing
end

function constructBC!(Ψ, n, D, Minv, bcleft, bcright, bctop, bcbottom, h1, h2, k1, k2)
    @inbounds @views Ψ1 = Ψ[3:(n - 1), 3:(n - 1)]
    @inbounds @views Ψ2 = Ψ[3:(n - 1), 2:n]'
    @inbounds @views Dcs = D[1, 3:(n - 1)]
    @inbounds @views Dce = D[n + 1, 3:(n - 1)]

    mul!(h1, Ψ1, Dcs)
    mul!(h2, Ψ1, Dce)

    @inbounds @views @. Ψ[3:(n - 1), 2] = Minv[1, 1] * (bctop[3:(n - 1)] - h1) +
                                          Minv[1, 2] * (bcbottom[3:(n - 1)] - h2)
    @inbounds @views @. Ψ[3:(n - 1), n] = Minv[2, 1] * (bctop[3:(n - 1)] - h1) +
                                          Minv[2, 2] * (bcbottom[3:(n - 1)] - h2)

    mul!(k1, Ψ2, Dcs)
    mul!(k2, Ψ2, Dce)

    @inbounds @views @. Ψ[2, 2:n] = Minv[1, 1] * (bcright[2:n] - k1) +
                                    Minv[1, 2] * (bcleft[2:n] - k2)
    @inbounds @views @. Ψ[n, 2:n] = Minv[2, 1] * (bcright[2:n] - k1) +
                                    Minv[2, 2] * (bcleft[2:n] - k2)

    return nothing
end

function f!(fu, u, p::CavityCache)
    @unpack fΨ,
    Ψ,
    n,
    Re,
    D1,
    D2,
    D4,
    bcleft,
    bcright,
    bctop,
    bcbottom,
    Minv,
    h1,
    h2,
    k1,
    k2,
    DΨ,
    ΨD,
    lDΨ,
    lΨD,
    laplΨ,
    biharmΨ,
    nonlinΨ = p

    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u

    constructΨ!(Ψ, n, D1, Minv, bcleft, bcright, bctop, bcbottom, h1, h2, k1, k2)

    mul!(DΨ, D4, Ψ)
    mul!(ΨD, Ψ, D4')

    mul!(laplΨ, D2, Ψ) # use as intermediate memory
    mul!(nonlinΨ, laplΨ, D2') # use as intermediate memory
    @. biharmΨ = DΨ + ΨD + 2 * nonlinΨ

    mul!(DΨ, D2, Ψ)
    mul!(ΨD, Ψ, D2')
    @. laplΨ = DΨ + ΨD

    mul!(DΨ, D1, Ψ)
    mul!(ΨD, Ψ, D1')

    mul!(lΨD, laplΨ, D1')
    mul!(lDΨ, D1, laplΨ)

    @. nonlinΨ = DΨ * lΨD - lDΨ * ΨD
    @. fΨ = (1 / Re) * biharmΨ - nonlinΨ
    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]

    return nothing
end

function f(u, p::CavityCache)
    dim = size(u, 1)
    fu = zeros(dim)
    f!(fu, u, p)

    return fu
end

function newton_custom!(f!, x0; tolmax = 1e-10, maxiter = 100)
    x = copy(x0)
    dim = size(x, 1)

    fx = similar(x)
    f!(fx, x)

    J = zeros(dim, dim)
    dx = zeros(dim)
    cache = FiniteDiff.JacobianCache(x)

    iter = 0
    tol = 1.0
    while tol > tolmax && iter < maxiter
        @time FiniteDiff.finite_difference_jacobian!(J, f!, x, cache)
        dx .= J \ (-fx)
        @. x = x + dx

        f!(fx, x)

        tol = norm(dx)
        iter += 1
    end

    return x, iter, tol
end
