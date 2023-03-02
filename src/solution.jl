using RecipesBase

abstract type Solution end

struct Solution1D <: Solution 
    nodes::Vector
    vals::Vector
end

RecipesBase.@recipe function f(sol::Solution1D)
    xlabel --> "nodes"
    ylabel --> "u"
    return sol.nodes, sol.vals
end
