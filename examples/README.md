# Examples Cavity Flow

The example folder provides its own `Project.toml` file and the
`NS2DBenchmarkSolver` module is accessed as external package to illustrate the
usage.

The package has to be added locally using a relative path. In the Pkg REPL type :
```julia
pkg>activate .
pkg>add "../"
pkg>instantiate 
```

Checkout the package [Revise.jl](https://github.com/timholy/Revise.jl.git) if
you want to make changes to the package and then apply the changes to the
examples.

3 simple examples for the 4 sided cavity flow are listed below :
| --- | --- |
| example\_solve.jl | simple solve for steady, unstable, symmectric solution |
| example\_timestepping.jl | time-stepping starting with random noised (saved as a gif)  |
| example\_continuation.jl | Reynolds continuation |
