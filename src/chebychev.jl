"""
    chebychev(N)

Compute the one-dimensional Chebichev Grid `x` using `N` points and 
the associated Chebychev Differentation matrix.

# Example
```jldoctest
julia> chebychev(2)
[1.5 -2.0 0.5; 0.5 0.0 -0.5; -0.5 2.0 -1.5], [1.0; 0.0; -1.0]
```
"""
function chebychev(N::Int)
    if N == 0
        D = [0]
        x = [1.0]
        return D, x
    end
    x = cos.(pi*(0:N)/N)
    c = [2; ones(N-1); 2].*(-1).^(0:N)
    X = repeat(x, 1, N+1)
    dX = X .- X'
    D = (c * (1.0./c)') ./ (dX + I)
    D = D - diagm(vec(sum(D, dims=2)))
    return D, x
end
