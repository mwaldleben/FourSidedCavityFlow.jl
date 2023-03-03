using RecipesBase

RecipesBase.@recipe function f(sol::Solution1D)
    xlabel --> "nodes"
    ylabel --> "u"
    return sol.nodes, sol.vals
end
