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

n = 64
p = CF.setup_struct(n, 0)

println("--- Run $(n)x$(n): ---")

# Folder where results are stored 
folder = "$(n)x$(n)"
mkdir(folder)
include("helpers.jl")

# Continuation
foldercont = "$folder/cont"
mkdir(foldercont)
include("study_continuation.jl")

# Linear stability analysis
folderlsa = "$folder/lsa"
mkdir(folderlsa)
include("study_linearstability.jl")

# Periodic orbits
folderpo = "$folder/po"
mkdir(folderpo)
include("study_periodicorbits.jl")

# Second branch
folderbranch2 = "$folder/branch2"
folderswitch = "$folderbranch2/switch_branch"
foldercont_branch2 = "$folderbranch2/cont"
mkdir(folderbranch2)
mkdir(folderswitch)
mkdir(foldercont_branch2)
include("study_branch2.jl")
