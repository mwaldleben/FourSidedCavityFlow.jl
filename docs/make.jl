push!(LOAD_PATH, "../")
push!(LOAD_PATH, "../src")
using Documenter, FourSidedCavityFlow

makedocs(;
    doctest = false,
    sitename = "FourSidedCavityFlow.jl",
    modules = [FourSidedCavityFlow],
    pages = ["Home" => "index.md",
                "Four-Sided Cavity Flows" => "foursidedcavityflows.md",
                "Tutorial" => "tutorial.md",
                "Reproduce Results" => "reproducibility.md",
                "API of Functions" => "api.md",
            ],
    format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
        )
    )

deploydocs(; repo = "github.com/morwald/FourSidedCavityFlow.jl.git")
