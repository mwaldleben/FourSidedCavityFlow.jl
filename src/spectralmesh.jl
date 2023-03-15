function SpectralMesh1D(n::Int, length::Real)
    nodes, D1, D2, D4  = diffchebyshev(n, length=length)

    return SpectralMesh1D(n, length, nodes, D1, D2, D4)
end

function SpectralMesh1D(n::Int)
    length = 1

    nodes, D1, D2, D4  = diffchebyshev(n)

    return SpectralMesh1D(n, length, nodes, D1, D2, D4)
end

function SpectralMesh2D(nx::Int, ny::Int, xlength::Real, ylength::Real)
    xnodes, Dx1, Dx2, Dx4  = diffchebyshev(nx, length=xlength)
    ynodes, Dy1, Dy2, Dy4  = diffchebyshev(ny, length=xlength)

    return SpectralMesh2D(nx, ny, xlength, ylength, xnodes, ynodes, Dx1, Dx2, Dx4, Dy1, Dy2, Dy4)
end

function SpectralMesh2D(n::Tuple{Int, Int})
    xlength = 1
    ylength = 1

    xnodes, Dx1, Dx2, Dx4  = diffchebyshev(n[1])
    ynodes, Dy1, Dy2, Dy4  = diffchebyshev(n[2])

    return SpectralMesh2D(n[1], n[2], xlength, ylength, xnodes, ynodes, Dx1, Dx2, Dx4, Dy1, Dy2, Dy4)
end
