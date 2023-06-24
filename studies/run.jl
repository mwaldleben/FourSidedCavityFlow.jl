using FourSidedCavityFlow
using DelimitedFiles
using CSV
using DataFrames
using Printf
using UnPack
using LaTeXStrings
using Random
const CF = FourSidedCavityFlow

# Random.seed!(1234)

n = 32
p = CF.setup_struct(n, 0)

println("--- Run $(n)x$(n): ---")

# Folder where results are stored 
folder = "$(n)x$(n)"
mkdir(folder)
include("helpers.jl")

foldercont = "$folder/cont"
mkdir(foldercont)
include("study_continuation.jl")

folderlsa = "$folder/lsa"
mkdir(folderlsa)
include("study_linearstability.jl")

# folderbranch2 = "$folder/branch2"
# folderswitch = "$folderbranch2/switch_branch"
# foldercont_branch2 = "$folderbranch2/cont"
# mkdir(folderbranch2)
# mkdir(folderswitch)
# mkdir(foldercont_branch2)
# include("study_branch2.jl")
