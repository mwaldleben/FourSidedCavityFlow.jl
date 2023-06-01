mutable struct CavityParameters{T <: Real}
    # Parameters which are change
    Re::T
    Ψstart::Matrix{T}

    # Normally fixed parameters 
    n::Int64
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
    h1
    h2
    k1
    k2
    fΨ
    Ψ
    D2Ψ
    ΨD2
    D4Ψ
    ΨD4
    biharmΨ
    laplΨ
    Ψ0
    laplΨ0
    nonlinΨ
end

struct CavityStruct{T <: Real}
    params::CavityParameters{T}
    cache::CavityCache
end

function setup_struct(n, Re)
    nodes, D1, D2, D4 = diff_chebyshev(n)

    k0 = 10
    bcfunc(x) = ((exp(k0 * (x - 1)) - 1) * (exp(-k0 * (x + 1)) - 1))^2

    bcleft = -bcfunc.(nodes)
    bcright = bcfunc.(nodes)
    bctop = bcfunc.(nodes)
    bcbottom = -bcfunc.(nodes)
    Minv = constructBC_matrix(D1)

    scl = 1e6
    Ψstart = zeros(n + 1, n + 1)

    params = CavityParameters{Float64}(Re, Ψstart, n, nodes, D1, D2, D4, bcleft, bcright, bctop,
                                       bcbottom, Minv[1, 1], Minv[1, 2], Minv[2, 1], Minv[2, 2], scl)

    # h1 = similar(bctop[3:(n - 1)])
    # h2 = similar(bctop[3:(n - 1)])
    # k1 = similar(bctop[2:n])
    # k2 = similar(bctop[2:n])
    #
    # Ψ = zeros(n + 1, n + 1)
    # fΨ = similar(Ψ)
    # D2Ψ = similar(Ψ)
    # ΨD2 = similar(Ψ)
    # D4Ψ = similar(Ψ)
    # ΨD4 = similar(Ψ)
    # laplΨ = similar(Ψ)
    # biharmΨ = similar(Ψ)
    # Ψ0 = zeros(n + 1, n + 1)
    # laplΨ0 = similar(Ψ)
    # nonlinΨ = similar(Ψ)

    vectmp = zeros(n - 3)
    vectmp2 = zeros(n - 3)
    vectmp3 = zeros(n - 1)
    vectmp4 = zeros(n - 1)
    tmp = zeros(n + 1, n + 1)
    tmp1 = zeros(n + 1, n + 1)
    tmp2 = zeros(n + 1, n + 1)
    tmp3 = zeros(n + 1, n + 1)
    tmp4 = zeros(n + 1, n + 1)
    tmp5 = zeros(n + 1, n + 1)
    tmp6 = zeros(n + 1, n + 1)
    tmp7 = zeros(n + 1, n + 1)
    tmp8 = zeros(n + 1, n + 1)
    tmp9 = zeros(n + 1, n + 1)
    tmp10 = zeros(n + 1, n + 1)

    h1 = DiffCache(vectmp, 12)
    h2 = DiffCache(vectmp2, 12)
    k1 = DiffCache(vectmp3, 12)
    k2 = DiffCache(vectmp4, 12)
    fΨ = DiffCache(tmp, 12)
    Ψ = DiffCache(tmp1, 12)
    laplΨ = DiffCache(tmp2, 12)
    biharmΨ = DiffCache(tmp3, 12)
    D2Ψ = DiffCache(tmp4, 12)
    ΨD2 = DiffCache(tmp5, 12)
    D4Ψ = DiffCache(tmp6, 12)
    ΨD4 = DiffCache(tmp7, 12)
    Ψ0 = DiffCache(tmp8, 12)
    laplΨ0 = DiffCache(tmp9, 12)
    nonlinΨ = DiffCache(tmp10, 12)

    cache = CavityCache{Float64}(h1, h2, k1, k2, Ψ, fΨ, D2Ψ, ΨD2, D4Ψ, ΨD4,
                                 biharmΨ, laplΨ, Ψ0, laplΨ0, nonlinΨ)

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

    h1 = get_tmp(h1, Ψ)
    h2 = get_tmp(h2, Ψ)
    k1 = get_tmp(k1, Ψ)
    k2 = get_tmp(k2, Ψ)

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

    h1 = get_tmp(h1, Ψ)
    h2 = get_tmp(h2, Ψ)
    k1 = get_tmp(k1, Ψ)
    k2 = get_tmp(k2, Ψ)

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

    Ψ = get_tmp(Ψ, u)
    fΨ = get_tmp(fΨ, u)
    D2Ψ = get_tmp(D2Ψ, u)
    ΨD2 = get_tmp(ΨD2, u)
    D4Ψ = get_tmp(D4Ψ, u)
    ΨD4 = get_tmp(ΨD4, u)
    laplΨ = get_tmp(laplΨ, u)
    biharmΨ = get_tmp(biharmΨ, u)

    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u
    # @inbounds Ψ[3:(n - 1),3:(n - 1)] = reshape(u, n - 3 , n - 3)

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
    @inbounds @. fΨ = (1 / Re) * biharmΨ - laplΨ

    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]
    # @inbounds fu .= reshape(fΨ[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    return nothing
end

function ftime!(fu, u, p::CavityStruct, Δt)
    @unpack Re, Ψstart, n, D1, D2, D4 = p.params
    @unpack fΨ, Ψ, D2Ψ, ΨD2, D4Ψ, ΨD4, laplΨ, biharmΨ, laplΨ0, nonlinΨ = p.cache

    Ψ = get_tmp(Ψ, u)
    fΨ = get_tmp(fΨ, u)
    D2Ψ = get_tmp(D2Ψ, u)
    ΨD2 = get_tmp(ΨD2, u)
    D4Ψ = get_tmp(D4Ψ, u)
    ΨD4 = get_tmp(ΨD4, u)
    laplΨ = get_tmp(laplΨ, u)
    biharmΨ = get_tmp(biharmΨ, u)
    laplΨ0 = get_tmp(laplΨ0, u)
    nonlinΨ = get_tmp(nonlinΨ, u)

    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u
    # @inbounds Ψ[3:(n - 1), 3:(n - 1)] = reshape(u, n - 3, n - 3)

    constructBC!(Ψ, p)

    mul!(D4Ψ, D4, Ψ)
    mul!(ΨD4, Ψ, D4')

    mul!(D2Ψ, D2, Ψ)
    mul!(ΨD2, Ψ, D2')

    mul!(laplΨ, D2Ψ, D2') # using as intermediate memory

    @inbounds @. biharmΨ = D4Ψ + ΨD4 + 2 * laplΨ
    @inbounds @. laplΨ = D2Ψ + ΨD2

    mul!(D2Ψ, D2, Ψstart)
    mul!(ΨD2, Ψstart, D2')

    @inbounds @. laplΨ0 = D2Ψ + ΨD2

    mul!(D2Ψ, D1, Ψ)
    mul!(ΨD2, Ψ, D1')

    mul!(ΨD4, laplΨ, D1')
    mul!(D4Ψ, D1, laplΨ)

    @inbounds @. nonlinΨ = D2Ψ * ΨD4 - D4Ψ * ΨD2
    @inbounds @. fΨ = (1 / Re) * biharmΨ - nonlinΨ

    @inbounds @. fΨ = Δt * fΨ - laplΨ + laplΨ0

    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]
    # @inbounds fu .= reshape(fΨ[3:(n - 1), 3:(n - 1)], (n - 3) * (n - 3))

    return nothing
end
