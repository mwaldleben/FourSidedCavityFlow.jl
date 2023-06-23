mutable struct CavityParameters{T <: Real}
    # Parameters which are change
    Re::T
    Ψi::Matrix{T}

    # Fixed parameters 
    n::Int64  # dimension: n+1 
    ic::Int64 # index center value
    i1::Int64 # index quarter 
    i2::Int64 # index three quarters
    nodes::Vector{T}
    D1::Matrix{T}
    D2::Matrix{T}
    D4::Matrix{T}
    bcleft::Vector{T}
    bcright::Vector{T}
    bctop::Vector{T}
    bcbottom::Vector{T}
    m11::T
    m12::T
    m21::T
    m22::T
    scl::T
end

struct CavityCache{T <: Real}
    h1::Vector{T}
    h2::Vector{T}
    k1::Vector{T}
    k2::Vector{T}
    Ψ::Matrix{T}
    fΨ::Matrix{T}
    D2Ψ::Matrix{T}
    ΨD2::Matrix{T}
    D4Ψ::Matrix{T}
    ΨD4::Matrix{T}
    biharmΨ::Matrix{T}
    laplΨ::Matrix{T}
    Ψ0::Matrix{T}
    laplΨ0::Matrix{T}
    nonlinΨ::Matrix{T}
end

struct CavityStruct{T <: Real}
    params::CavityParameters{T}
    cache::CavityCache{T}
end

function setup_struct(n, Re)
    nodes, D1, D2, D4 = diff_chebyshev(n)

    ic = floor(n / 2) + 1
    i1 = floor(n / 4) + 1
    i2 = 3 * floor(n / 4) + 1

    k0 = 10
    bcfunc(x) = ((exp(k0 * (x - 1)) - 1) * (exp(-k0 * (x + 1)) - 1))^2

    bcleft = -bcfunc.(nodes)
    bcright = bcfunc.(nodes)
    bctop = bcfunc.(nodes)
    bcbottom = -bcfunc.(nodes)

    Minv = constructBC_matrix(D1)

    scl = 1e6
    Ψi = zeros(n + 1, n + 1)

    params = CavityParameters{Float64}(Re,
        Ψi,
        n,
        ic,
        i1,
        i2,
        nodes,
        D1,
        D2,
        D4,
        bcleft,
        bcright,
        bctop,
        bcbottom,
        Minv[1, 1],
        Minv[1, 2],
        Minv[2, 1],
        Minv[2, 2],
        scl)

    h1 = similar(bctop[3:(n - 1)])
    h2 = similar(bctop[3:(n - 1)])
    k1 = similar(bctop[2:n])
    k2 = similar(bctop[2:n])

    Ψ = zeros(n + 1, n + 1)
    fΨ = similar(Ψ)
    D2Ψ = similar(Ψ)
    ΨD2 = similar(Ψ)
    D4Ψ = similar(Ψ)
    ΨD4 = similar(Ψ)
    laplΨ = similar(Ψ)
    biharmΨ = similar(Ψ)
    Ψ0 = zeros(n + 1, n + 1)
    laplΨ0 = similar(Ψ)
    nonlinΨ = similar(Ψ)

    cache = CavityCache{Float64}(h1, h2, k1, k2, Ψ, fΨ, D2Ψ, ΨD2, D4Ψ, ΨD4, biharmΨ, laplΨ,
        Ψ0, laplΨ0, nonlinΨ)

    p = CavityStruct{Float64}(params, cache)

    return p
end

function constructBC_matrix(D)
    n = size(D, 1) - 1
    M = [D[1, 2] D[1, n]; D[n + 1, 2] D[n + 1, n]]
    return inv(M)
end

function constructBC!(Ψ, p::CavityStruct)
    @unpack n, D1, m11, m12, m21, m22, bcleft, bcright, bctop, bcbottom = p.params
    @unpack h1, h2, k1, k2 = p.cache

    @inbounds @views Ψ1 = Ψ[3:(n - 1), 3:(n - 1)]
    @inbounds @views Ψ2 = Ψ[3:(n - 1), 2:n]'
    @inbounds @views Dcs = D1[1, 3:(n - 1)]
    @inbounds @views Dce = D1[n + 1, 3:(n - 1)]

    mul!(h1, Ψ1, Dcs)
    mul!(h2, Ψ1, Dce)

    @inbounds @views @. Ψ[3:(n - 1), 2] = m11 * (bctop[3:(n - 1)] - h1) +
                                          m12 * (bcbottom[3:(n - 1)] - h2)
    @inbounds @views @. Ψ[3:(n - 1), n] = m21 * (bctop[3:(n - 1)] - h1) +
                                          m22 * (bcbottom[3:(n - 1)] - h2)

    mul!(k1, Ψ2, Dcs)
    mul!(k2, Ψ2, Dce)

    @inbounds @views @. Ψ[2, 2:n] = m11 * (bcleft[2:n] - k1) + m12 * (bcright[2:n] - k2)
    @inbounds @views @. Ψ[n, 2:n] = m21 * (bcleft[2:n] - k1) + m22 * (bcright[2:n] - k2)

    return nothing
end

function constructBC(u, p::CavityStruct)
    @unpack n = p.params

    Ψ = zeros(n + 1, n + 1)
    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u

    constructBC!(Ψ, p)

    return Ψ
end

function construct_homogenousBC!(Ψ, p::CavityStruct)
    @unpack n, D1, m11, m12, m21, m22, bcleft, bcright, bctop, bcbottom = p.params
    @unpack h1, h2, k1, k2 = p.cache

    @inbounds @views Ψ1 = Ψ[3:(n - 1), 3:(n - 1)]
    @inbounds @views Ψ2 = Ψ[3:(n - 1), 2:n]'
    @inbounds @views Dcs = D1[1, 3:(n - 1)]
    @inbounds @views Dce = D1[n + 1, 3:(n - 1)]

    mul!(h1, Ψ1, Dcs)
    mul!(h2, Ψ1, Dce)

    @inbounds @views @. Ψ[3:(n - 1), 2] = -m11 * h1 - m12 * h2
    @inbounds @views @. Ψ[3:(n - 1), n] = -m21 * h1 - m22 * h2

    mul!(k1, Ψ2, Dcs)
    mul!(k2, Ψ2, Dce)

    @inbounds @views @. Ψ[2, 2:n] = -m11 * k1 - m12 * k2
    @inbounds @views @. Ψ[n, 2:n] = -m21 * k1 - m22 * k2
    return nothing
end

function f!(fu, u, p::CavityStruct)
    @unpack Re, n, D1, D2, D4 = p.params
    @unpack fΨ, Ψ, D2Ψ, ΨD2, D4Ψ, ΨD4, laplΨ, biharmΨ = p.cache

    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u

    constructBC!(Ψ, p)

    mul!(D4Ψ, D4, Ψ)
    mul!(ΨD4, Ψ, D4')

    mul!(D2Ψ, D2, Ψ)
    mul!(ΨD2, Ψ, D2')

    mul!(laplΨ, D2Ψ, D2') # using as intermediate memory
    @inbounds @. biharmΨ = D4Ψ + ΨD4 + 2 * laplΨ
    @inbounds @. laplΨ = D2Ψ + ΨD2

    mul!(D2Ψ, D1, Ψ)
    mul!(ΨD2, Ψ, D1')

    mul!(ΨD4, laplΨ, D1')
    mul!(D4Ψ, D1, laplΨ)

    @inbounds @. laplΨ = D2Ψ * ΨD4 - D4Ψ * ΨD2
    @inbounds @. fΨ = (1 / Re) * biharmΨ + laplΨ

    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]

    return nothing
end

function ftime!(fu, u, p::CavityStruct, Δt)
    @unpack Re, Ψi, n, D1, D2, D4 = p.params
    @unpack fΨ, Ψ, D2Ψ, ΨD2, D4Ψ, ΨD4, laplΨ, biharmΨ, laplΨ0, nonlinΨ = p.cache

    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u

    constructBC!(Ψ, p)

    mul!(D4Ψ, D4, Ψ)
    mul!(ΨD4, Ψ, D4')

    mul!(D2Ψ, D2, Ψ)
    mul!(ΨD2, Ψ, D2')

    mul!(laplΨ, D2Ψ, D2') # using as intermediate memory

    @inbounds @. biharmΨ = D4Ψ + ΨD4 + 2 * laplΨ
    @inbounds @. laplΨ = D2Ψ + ΨD2

    mul!(D2Ψ, D2, Ψi)
    mul!(ΨD2, Ψi, D2')

    @inbounds @. laplΨ0 = D2Ψ + ΨD2

    mul!(D2Ψ, D1, Ψ)
    mul!(ΨD2, Ψ, D1')

    mul!(ΨD4, laplΨ, D1')
    mul!(D4Ψ, D1, laplΨ)

    @inbounds @. nonlinΨ = D2Ψ * ΨD4 - D4Ψ * ΨD2
    @inbounds @. fΨ = (1 / Re) * biharmΨ + nonlinΨ

    @inbounds @. fΨ = Δt * fΨ - laplΨ + laplΨ0

    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]

    return nothing
end
