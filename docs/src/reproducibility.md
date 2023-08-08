# Reproduce Results

The `studies` folder in the project directory provides the necessary files to reproduce the results.
This folder provides its own `Project.toml` file and the `FourSidedCavityFlow`
module is accessed as an external package.

The module can be added locally using a relative path. In the `Pkg` REPL type :
```julia
pkg>activate .
pkg>dev "../"
pkg>instantiate 
```

| File                            | Description                                                 |
|:------------------------------- |:----------------------------------------------------------- |
| run.jl                          | Entry point, set the grid size, runs files listed below     |
| study_continuation.jl           | Run pseudo-arclength continuation                           |
| study_linearstability.jl        | Linear stability analysis on points of continuation         |
| study_branch2.jl                | Continuation of second branch and linear stability analysis |
| study_periodicorbits.jl         | Generate periodic orbits of the Hopf bifurcation            |
