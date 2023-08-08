# Studies for the Regularized Four-sided Cavity

This folder provides its own `Project.toml` file and the `FourSidedCavityFlow`
module is accessed as an external package to reproduce the results.

The module cab be added locally using a relative path. In the `Pkg` REPL type :
```julia
pkg>activate .
pkg>dev "../"
pkg>instantiate 
```

Check out the package [Revise.jl](https://github.com/timholy/Revise.jl.git) if
you want to make changes to the package and then apply the changes to the
studies.

| File                            | Description                                                 |
|:------------------------------- |:----------------------------------------------------------- |
| run.jl                          | Entry point, set the grid size, runs files listed below     |
| study_continuation.jl           | Run pseudo-arclength continuation                           |
| study_linearstability.jl        | Linear stability analysis on points of continuation         |
| study_branch2.jl                | Continuation of second branch and linear stability analysis |
| study_periodicorbits.jl         | Generate periodic orbits of the Hopf bifurcation            |
