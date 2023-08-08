"""
A `struct` containing the parameters for the four-sided cavity flow. The fixed
parameters consist of the fields for the discretization and the boundary conditions. The 
field `Re` is meant to be set beforehand to run a simulation for a given Reynolds number. 
"""
mutable struct CavityParameters{T <: Real}
    # Parameters which are change
    Re::T # Reynolds number
    Ψi::Matrix{T} # solution of previous timestep, when doing time-integration

    # Fixed parameters 
    n::Int64  # dimension: n+1 
    ic::Int64 # index center value
    i1::Int64 # index quarter 
    i2::Int64 # index three quarters
    nodes::Vector{T} # nodes in [1,-1] 
    D1::Matrix{T} # First Chebyshev differentation  
    D2::Matrix{T}
    D4::Matrix{T}
    bcleft::Vector{T} # streamfunction boundary conditions 
    bcright::Vector{T}
    bctop::Vector{T}
    bcbottom::Vector{T}
    m11::T # elements of 2x2 matrix for BC calculations
    m12::T
    m21::T
    m22::T
    scl::T # scaling for augmented system in pseudo-arclength continuation
end

"""
A `struct` containing the cache variables for the nonlinear function `f`.
"""
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

"""
The CavityStruct contains the sub-`struct`s CavityParameters and CavityCache. It is passed  
to to the function `f!(fu, u, p::CavityStruct)` to provide the parameters and the cache for the 
evaluation of the nonlinear function.
"""
struct CavityStruct{T <: Real}
    params::CavityParameters{T}
    cache::CavityCache{T}
end

@doc raw"""
    setup_struct(n, Re)

Function to instantiate the CavityStruct. Creates a grid of size ``(n+1) \times (n+1)``. The regularization 
parameter for the boundary conditions is fixed to ``k_0 = 10``.

# Example
```julia-repl
julia> p = setup_struct(64, 1)
```
"""
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

"""
    Minv = constructBC_matrix(D)

Helper function to return the 2 by 2 matrices needed for the construction of the boundary for the streamfunction.
`D` is the Chebyshev differentiation matrix of the discretization. 
"""
function constructBC_matrix(D)
    n = size(D, 1) - 1
    M = [D[1, 2] D[1, n]; D[n + 1, 2] D[n + 1, n]]
    return inv(M)
end

"""
    construct!BC(Ψ, p)

Constructs the two outer grid points rows and columns from the given boundary conditions. The outer rows and columns of
the `Ψ` matrix will be overridden.
"""
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

"""
    Ψ = constructBC(u, p)

Out-of-place version of the function constructBC!(Ψ, p). The input is the vector of inner grid points `u`
and the function returns the `Ψ` matrix with the boundaries reconstructed for the given flow field.
"""
function constructBC(u, p::CavityStruct)
    @unpack n = p.params

    Ψ = zeros(n + 1, n + 1)
    @inbounds @views Ψ[3:(n - 1), 3:(n - 1)][:] .= u

    constructBC!(Ψ, p)

    return Ψ
end

"""
    construct_homogenousBC(Ψ, p)

Constructs the boundaries for the streamfunction matrix `Ψ` by using homogenous boundary conditions
(zero velocities at all lids of the cavity). This is employed in the perturbed solution for the linear stability analysis.
"""
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

@doc raw"""
    f!(fu, u, p)

The nonlinear function for steady-state solutions of the four-sided cavity flow. Evaluates the equation:

``F(\Psi) = \frac{1}{\mathrm{Re}} \Delta^2 \Psi + (\partial_x \Psi) \partial_y(\Delta \Psi) - (\partial_y \Psi) \partial_x(\Delta \Psi)``

Only the interior part `u` (given as a vector) of the streamfunction matrix `Ψ` is needed, 
as the boundary conditions explicitly give the 2 outer rows and columns. `fu` is a vector where the solution of the function evaluation will be stored.

# Example
```julia-repl
julia> n = 64; Re = 100; 
julia> p = setup_struct(n, Re) 
julia> Ψ = zeros(n + 1, n+ 1)
julia> u = Ψ[3:(n - 1), 3:(n - 1)][:]
julia> fu = similiar(u) 
julia> f!(fu, u, p) 
```
"""
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


@doc raw"""
    ftime!(fu, u, p, h)

The nonlinear function for the time integration of the streamfunction formulation. Solves the equations discretized with an
implicit Euler scheme and a timestep `h`:

``F(\Psi_{i+1}) = \Delta\Psi_i - \Delta\Psi_{i+1} + h (\frac{1}{\mathrm{Re}}\Delta^2 \Psi_{i+1} +
  (\partial_x \Psi_{i+1}) \partial_y(\Delta \Psi_{i+1}) -
  (\partial_y \Psi_{i+1}) \partial_x(\Delta \Psi_{i+1})) = 0``

The solution of the previous timestep has to be set beforehand by `p.params.Ψi`.
"""
function ftime!(fu, u, p::CavityStruct, h)
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

    @inbounds @. fΨ = h * fΨ - laplΨ + laplΨ0

    @inbounds @views fu .= fΨ[3:(n - 1), 3:(n - 1)][:]

    return nothing
end
