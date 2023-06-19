"""
    diff_chebyshev(n, length)

Compute the Chebyshev differentiation matrix with `n`+1 points in the
interval [length, -length] which defaults to [1, -1]. Returns a vector of nodes 
and the matrices for the first, second and fourth derivatives.
"""
function diff_chebyshev(n::Int; length = 1::Real)
    # Define nodes in [1, -1]
    nodes = [cos(k * π / n) for k in 0:n]

    # Off-diagonal entries
    c = [2; ones(n - 1); 2]
    dij = (i, j) -> (-1)^(i + j) * c[i + 1] / (c[j + 1] * (nodes[i + 1] - nodes[j + 1]))
    D1 = [dij(i, j) for i in 0:n, j in 0:n]

    # Diagonal entries
    D1[isinf.(D1)] .= 0
    s = sum(D1; dims = 2)
    D1 -= diagm(s[:, 1])

    # Change interval if necessary
    if length ≠ 1
        a = length / 2
        b = -length / 2
        nodes = @. a + (b - a) * (nodes - 1) / 2
        D1 = 2 * D1 / (b - a)
    end

    # Second derivative
    D2 = D1^2

    # Fourth derivative
    D4 = D2^2

    return nodes, D1, D2, D4
end
