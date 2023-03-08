using RecipesBase

function print(sol::Solution2D)
    xnbcells = size(sol.xnodes, 1)
    ynbcells = size(sol.ynodes, 1)
    isconverged = sol.isconverged
    iter = sol.iter

    println("-----------------------------")
    println("mesh: $xnbcells x $ynbcells")
    println("converged: $isconverged")
    println("iterations: $iter")
    println("-----------------------------")
end

RecipesBase.@recipe function f(sol::Solution1D)
    xlabel --> "nodes"
    ylabel --> "u"
    return sol.nodes, sol.vals
end

RecipesBase.@recipe function f(sol::Solution2D)
    xlabel --> "x"
    ylabel --> "y"
    return sol.xnodes, sol.ynodes, sol.vals
end
