# Studies for the regularized four-sided cavity

This folder provides its own `Project.toml` file and the `FourSidedCavityFlow`
module is accessed as an external package to reproduce the results.

The package has to be added locally using a relative path. In the Pkg REPL type :
```julia
pkg>activate .
pkg>dev "../"
pkg>instantiate 
```

Check out the package [Revise.jl](https://github.com/timholy/Revise.jl.git) if
you want to make changes to the package and then apply the changes to the
examples.

| File                            | Description                                                |
| ------------------------------- | ---------------------------------------------------------- |
