"""
    diffchebychev(n, span)

Compute the Chebyshev differentiation matrix with `n`+1 points in the
interval `span` which defaults to [-1. 1]. Returns a vector of nodes 
and the matrices for the first, second and fourth derivatives.
"""
function diffchebychev(n::Int, span=[-1, 1])
    # Define nodes in [-1,1]
    x = [ -cos(k*π/n) for k in 0:n ]    

    # Off-diagonal entries
    c = [2; ones(n-1); 2];
    dij = (i,j) -> (-1)^(i+j)*c[i+1]/(c[j+1]*(x[i+1]-x[j+1]))
    Dₓ = [ dij(i,j) for i in 0:n, j in 0:n ]

    # Diagonal entries
    Dₓ[isinf.(Dₓ)] .= 0
    s = sum(Dₓ,dims=2)
    Dₓ -= diagm(s[:,1])

    # Change interval if necessary
    if span[1] ≠ -1 && span[2] ≠ 1 
        a,b = span
        x = @. a + (b-a)*(x+1)/2
        Dₓ = 2*Dₓ/(b-a)
    end

    # Second derivative
    Dₓₓ = Dₓ^2

    # Fourth derivative
    Dₓₓₓₓ = Dₓₓ^2

    return x, Dₓ, Dₓₓ, Dₓₓₓₓ
end
