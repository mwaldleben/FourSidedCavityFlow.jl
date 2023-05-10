mutable struct CavityParameters{T <: Real}
    # Parameters
    Re::T

    # Mesh and boundary paramaters 
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

    # Cache for calculations
    h1::Vector{T}
    h2::Vector{T}
    k1::Vector{T}
    k2::Vector{T}
    scl::T
    fΨ::Matrix{T}
    Ψ::Matrix{T}
    Ψ0::Matrix{T}
    D2Ψ::Matrix{T}
    ΨD2::Matrix{T}
    D4Ψ::Matrix{T}
    ΨD4::Matrix{T}
    biharmΨ::Matrix{T}
    laplΨ::Matrix{T}
    laplΨ0::Matrix{T}
    nonlinΨ::Matrix{T}
end

function construct_BC_matrix(D)
    n = size(D, 1) - 1
    M = [D[1, 2] D[1, n]; D[n + 1, 2] D[n + 1, n]]
    return inv(M)
end

function construct_BC!(p)
    @unpack n, Ψ, D1, m11, m12, m21, m22, bcleft, bcright, bctop, bcbottom, h1, h2, k1, k2 = p

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

function construct_BC(p)
    construct_BC!(p)

    return copy(p.Ψ)
end

function construct_homogenous_BC!(p)
    @unpack n, Ψ, D1, m11, m12, m21, m22, bcleft, bcright, bctop, bcbottom, h1, h2, k1, k2 = p

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

function f!(fu, u, p::CavityParameters)
    @unpack Re, n, D1, D2, D4, fΨ, Ψ, D2Ψ, ΨD2, D4Ψ, ΨD4, laplΨ, biharmΨ = p

    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u
    # @inbounds Ψ[3:(n - 1),3:(n - 1)] = reshape(u, n - 3 , n - 3)

    construct_BC!(p)

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

function ftime!(fu, u, p, Δt)
    @unpack Re, n, D1, D2, D4, fΨ, Ψ, Ψ0, D2Ψ, ΨD2, D4Ψ, ΨD4, laplΨ, biharmΨ, laplΨ0, nonlinΨ = p

    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u
    # @inbounds Ψ[3:(n - 1), 3:(n - 1)] = reshape(u, n - 3, n - 3)

    construct_BC!(p)

    mul!(D4Ψ, D4, Ψ)
    mul!(ΨD4, Ψ, D4')

    mul!(D2Ψ, D2, Ψ)
    mul!(ΨD2, Ψ, D2')

    mul!(laplΨ, D2Ψ, D2') # using as intermediate memory

    @inbounds @. biharmΨ = D4Ψ + ΨD4 + 2 * laplΨ
    @inbounds @. laplΨ = D2Ψ + ΨD2

    mul!(D2Ψ, D2, Ψ0)
    mul!(ΨD2, Ψ0, D2')

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

function setup_params(n, Re)
    nodes, D1, D2, D4 = diff_chebyshev(n)

    k0 = 10
    bcfunc(x) = ((exp(k0 * (x - 1)) - 1) * (exp(-k0 * (x + 1)) - 1))^2

    bcleft = -bcfunc.(nodes)
    bcright = bcfunc.(nodes)
    bctop = bcfunc.(nodes)
    bcbottom = -bcfunc.(nodes)

    Minv = construct_BC_matrix(D1)
    h1 = similar(bctop[3:(n - 1)])
    h2 = similar(bctop[3:(n - 1)])
    k1 = similar(bctop[2:n])
    k2 = similar(bctop[2:n])
    scl = 1e6

    Ψ = zeros(n + 1, n + 1)
    Ψ0 = zeros(n + 1, n + 1)
    fΨ = similar(Ψ)
    laplΨ = similar(Ψ)
    biharmΨ = similar(Ψ)
    laplΨ0 = similar(Ψ)
    nonlinΨ = similar(Ψ)
    D2Ψ = similar(Ψ)
    ΨD2 = similar(Ψ)
    D4Ψ = similar(Ψ)
    ΨD4 = similar(Ψ)

    p = CavityParameters{Float64}(Re, n, nodes, D1, D2, D4, bcleft, bcright, bctop,
                                  bcbottom, Minv[1, 1], Minv[1, 2], Minv[2, 1], Minv[2, 2],
                                  h1, h2, k1, k2, scl, fΨ, Ψ, Ψ0, D2Ψ, ΨD2, D4Ψ, ΨD4,
                                  biharmΨ, laplΨ, laplΨ0, nonlinΨ)

    return p
end
