"""
    diffchebychev(n, span)

Compute the Chebyshev differentiation matrix with `n`+1 points in the
interval `span` which defaults to [-1. 1]. Returns a vector of nodes 
and the matrices for the first, second and fourth derivatives.
"""
function diffchebychev(n::Int, span=(-1, 1))
    # Define nodes in [-1,1]
    x = [ -cos(k*π/n) for k in 0:n ]    

    # Off-diagonal entries
    c = [2; ones(n-1); 2];
    dij = (i,j) -> (-1)^(i+j)*c[i+1]/(c[j+1]*(x[i+1]-x[j+1]))
    Dx1 = [ dij(i,j) for i in 0:n, j in 0:n ]

    # Diagonal entries
    Dx1[isinf.(Dx1)] .= 0
    s = sum(Dx1,dims=2)
    Dx1 -= diagm(s[:,1])

    # Change interval if necessary
    if span[1] ≠ -1 && span[2] ≠ 1 
        a,b = span
        x = @. a + (b-a)*(x+1)/2
        Dx1 = 2*Dx1/(b-a)
    end

    # Second derivative
    Dx2 = Dx1^2

    # Fourth derivative
    Dx4 = Dx2^2

    return x, Dx1, Dx2, Dx4
end
