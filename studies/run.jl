using FourSidedCavityFlow
using DelimitedFiles
using Plots
using CSV
using DataFrames
using Printf
using UnPack
using LaTeXStrings
using Random
const CF = FourSidedCavityFlow
gr()

Random.seed!(1234)

n = 32
p = CF.setup_struct(n, 0)

println("--- Run $(n)x$(n): ---")

# Folder where results are stored 
foldercont = "continuation$(n)x$(n)"
# mkdir(foldercont)

# include("study_continuation.jl")
# include("study_converge_psis.jl")
# include("study_bifurcation_diag.jl")
# include("study_branch2.jl")
include("study_linearstability.jl")
include("study_converge_psis_bifurcation.jl")
