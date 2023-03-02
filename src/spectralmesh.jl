abstract type SpectralMesh end

struct SpectralMesh1D <: SpectralMesh
    nbcells::Int
    coordmin::Real
    coordmax::Real
    nodes::Vector
    diff1mat::Matrix
    diff2mat::Matrix
    diff4mat::Matrix
    spectralmethod::String
end

function SpectralMesh1D(nbcells, (coordmin, coordmax), spectralmethod="chebychev")
    nodes, diff1mat, diff2mat, diff4mat  = diffchebychev(nbcells, (coordmin, coordmax))

    return SpectralMesh1D(nbcells, coordmin, coordmax, nodes, diff1mat, diff2mat, diff4mat, spectralmethod)
end
