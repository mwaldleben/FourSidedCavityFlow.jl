function SpectralMesh1D(nbcells, span, spectralmethod="chebychev")
    nodes, diff1mat, diff2mat, diff4mat  = diffchebychev(nbcells, span=span)

    return SpectralMesh1D(nbcells, span[1], span[2], nodes, diff1mat, diff2mat, diff4mat, spectralmethod)
end

function SpectralMesh2D(nbcells, xspan, yspan, spectralmethod="chebychev")
    xnodes, diffx1mat, diffx2mat, diffx4mat  = diffchebychev(nbcells[1], span=xspan)
    ynodes, diffy1mat, diffy2mat, diffy4mat  = diffchebychev(nbcells[2], span=xspan)

    return SpectralMesh2D(nbcells[1], nbcells[2], xspan[1], xspan[2], yspan[1], yspan[2], xnodes, ynodes, diffx1mat, diffx2mat, diffx4mat, diffy1mat, diffy2mat, diffy4mat, spectralmethod)
end

function SpectralMesh2D(nbcells)
    xnodes, diffx1mat, diffx2mat, diffx4mat  = diffchebychev(nbcells[1])
    ynodes, diffy1mat, diffy2mat, diffy4mat  = diffchebychev(nbcells[2])
    spectralmethod = "chebychev"

    return SpectralMesh2D(nbcells[1], nbcells[2], -1, 1, -1, 1, xnodes, ynodes, diffx1mat, diffx2mat, diffx4mat, diffy1mat, diffy2mat, diffy4mat, spectralmethod)
end
