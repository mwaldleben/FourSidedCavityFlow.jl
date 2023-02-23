push!(LOAD_PATH,"../src/")
using Documenter, NS2DBenchmarkSolver

DocMeta.setdocmeta!(NS2DBenchmarkSolver, :DocTestSetup, :(using NS2DBenchmarkSolver); recursive=true)

makedocs(
         sitename = "NS2DBenchmarkSolver.jl",
         modules  = [NS2DBenchmarkSolver],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(
    repo = "github.com/morwald/NS2DBenchmarkSolver.jl.git",
)
