# FourSidedCavityFlow.jl Documentation

```@meta
CurrentModule = FourSidedCavityFlow
```

The four-sided cavity flow is a two-dimensional flow problem. It is an extension of the
simple one-sided lid-driven case, where all lids move with the same velocity
profile and parallel lids move in opposite directions.

```@raw html
<center><img src="./assets/foursidedcavity.png"></center>
```

```math
\partial_t \Delta \Psi = \frac{1}{\mathrm{Re}} \Delta^2 \Psi
  + (\partial_x \Psi) \partial_y(\Delta \Psi)
  - (\partial_y \Psi) \partial_x(\Delta \Psi)
```

This Julia module explores a regularized version of the four-sided lid-driven
cavity for incompressible fluids to be used as a validator benchmark for
Navier-Stokes solvers. The regularization overcomes the corner singularities
which are due to the discontinuous boundary conditions. The considered method
recovers exponential convergence with a pseudo-spectral Chebyshev
discretization scheme. Below the regularization function is shown.  

## Installation

This module is not a registered package. To install the FourSidedCavityFlow.jl, 
run the following commands in your shell.
```bash
git clone https://github.com/morwald/FourSidedCavityFlow.jl.git
cd FourSidedCavityFlow.jl
julia
```

Inside the Julia REPL open the built-in package manager Pkg by pressing `]` and then
run.
```julia
pkg>activate .
pkg>instantiate
```
This will activate the package and download the necessary dependencies.
