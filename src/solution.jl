using RecipesBase

RecipesBase.@recipe function f(sol::Solution1D)
    xlabel --> "nodes"
    ylabel --> "u"
    return sol.nodes, sol.vals
end

RecipesBase.@recipe function f(sol::Solution2D)
    xlabel --> ""
    ylabel --> "u"
    return sol.xnodes, sol.ynodes, sol.vals
end
