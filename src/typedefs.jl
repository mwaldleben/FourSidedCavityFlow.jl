mutable struct CavityParameters{T<:Real}
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
    fΨ::Matrix{T}
    Ψ::Matrix{T}
    Ψ0::Matrix{T}
    D2Ψ::Matrix{T}
    ΨD2::Matrix{T}
    D4Ψ::Matrix{T}
    ΨD4::Matrix{T}
    biharmΨ::Matrix{T}
    laplΨ::Matrix{T}
end
