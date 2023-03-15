# Examples of the four-sided cavity

The example folder provides its own `Project.toml` file and the
`CavityFlow` module is accessed as an external package to illustrate the
usage.

The package has to be added locally using a relative path. In the Pkg REPL type :
```julia
pkg>activate .
pkg>add "../"
pkg>instantiate 
```

Check out the package [Revise.jl](https://github.com/timholy/Revise.jl.git) if
you want to make changes to the package and then apply the changes to the
examples.

3 simple examples of the 4 sided cavity flow are listed below :

| File                            | Description                                                |
| ------------------------------- | ---------------------------------------------------------- |
| example\_solve.jl               | simple solve for steady, unstable, symmetric solution     |
| example\_timestepping.jl        | time-stepping starting with random noised (saved as a gif) |
| example\_naturalcontinuation.jl | natural continuation                                       |
