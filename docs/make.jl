push!(LOAD_PATH,"../src/")
using Documenter, CavityFlow

DocMeta.setdocmeta!(CavityFlow, :DocTestSetup, :(using CavityFlow); recursive=true)

makedocs(
         sitename = "CavityFlow.jl",
         modules  = [CavityFlow],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(
    repo = "github.com/morwald/CavityFlow.jl.git",
)
