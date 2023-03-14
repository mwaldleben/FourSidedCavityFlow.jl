function SpectralMesh1D(n::Int, length::Real)
    nodes, diff1, diff2, diff4  = diffchebychev(n, length=length)

    return SpectralMesh1D(n, length, nodes, diff1, diff2, diff4)
end

function SpectralMesh1D(n::Int)
    length = 1

    nodes, diff1, diff2, diff4  = diffchebychev(n)

    return SpectralMesh1D(n, length, nodes, diff1, diff2, diff4)
end

function SpectralMesh2D(nx::Int, ny::Int, lengthx::Real, lengthy::Real)
    nodesx, diffx1, diffx2, diffx4  = diffchebychev(nx, length=lengthx)
    nodesy, diffy1, diffy2, diffy4  = diffchebychev(ny, length=lengthy)

    return SpectralMesh2D(nx, ny, lengthx, lengthy, nodesx, nodesy, diffx1, diffx2, diffx4, diffy1, diffy2, diffy4)
end

function SpectralMesh2D(n::Tuple{Int, Int})
    lengthx = 1
    lengthy = 1

    nodesx, diffx1, diffx2, diffx4  = diffchebychev(n[1])
    nodesy, diffy1, diffy2, diffy4  = diffchebychev(n[2])

    return SpectralMesh2D(n[1], n[2], lengthx, lengthy, nodesx, nodesy, diffx1, diffx2, diffx4, diffy1, diffy2, diffy4)
end
