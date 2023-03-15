function setNeumannBC2D(probl::Cavity4Sided, bcleft::Vector, bcright::Vector, bctop::Vector, bcbottom::Vector)
    probl.bcleft = bcleft
    probl.bcright = bcright
    probl.bctop = bctop
    probl.bcbottom = bcbottom
end

function setNeumannBC2D(probl::Cavity4Sided, bcleft::Function, bcright::Function, bctop::Function, bcbottom::Function)
    probl.bcright = bcleft(probl.mesh.ynodes)
    probl.bcleft = bcright(probl.mesh.ynodes)
    probl.bctop = bctop(probl.mesh.xnodes)
    probl.bcbottom = bcbottom(probl.mesh.xnodes)
end
