push!(LOAD_PATH, "../src/")
using Documenter, FourSidedCavityFlow

DocMeta.setdocmeta!(FourSidedCavityFlow, :DocTestSetup, :(using FourSidedCavityFlow); recursive = true)

makedocs(sitename = "FourSidedCavityFlow.jl", modules = [FourSidedCavityFlow], pages = ["Home" => "index.md"])

deploydocs(repo = "github.com/morwald/FourSidedCavityFlow.jl.git")
