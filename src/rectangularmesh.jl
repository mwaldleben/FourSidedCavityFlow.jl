function ChebyshevMesh(nx::Int, ny::Int, xlength::Real, ylength::Real)
    xnodes, Dx1, Dx2, Dx4  = diffchebyshev(nx, length=xlength)
    ynodes, Dy1, Dy2, Dy4  = diffchebyshev(ny, length=xlength)

    return ChebyshevMesh(nx, ny, xlength, ylength, xnodes, ynodes, Dx1, Dx2, Dx4, Dy1, Dy2, Dy4)
end

function ChebyshevMesh(n::Tuple{Int, Int})
    xlength = 1
    ylength = 1

    xnodes, Dx1, Dx2, Dx4  = diffchebyshev(n[1])
    ynodes, Dy1, Dy2, Dy4  = diffchebyshev(n[2])

    return ChebyshevMesh(n[1], n[2], xlength, ylength, xnodes, ynodes, Dx1, Dx2, Dx4, Dy1, Dy2, Dy4)
end
