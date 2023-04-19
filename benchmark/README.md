# Benchmark

Compare solve steady state with a precomputed solution from timestepping at Reynolds 100.
Median values of times are compared.

| Test                   | Grid    | Matlab ex. time (s) | Julia ex. time (s) |
| ---------------------- | ------- | ------------------- | ------------------ |
| Function evaluation    | 32 x 32 | 0.000094894         | 0.000052603        |  
|                        | 64 x 64 | 0.00032666          | 0.000223192        |
| Creating the jacobian  | 32 x 32 | 0.0836              | 0.049540           |
|                        | 64 x 64 | 1.1526              | 0.858283           |
| Newton (3 iterations)  | 32 x 32 | 0.1967              | 0.126673           |
|                        | 64 x 64 | 5.8606              | 4.878              |

# Run script from the command line
```bash
julia --project=CavityFlow.jl/benchmark CavityFlow.jl/benchmark/benchmark_solve.jl > benchmark_results.txt
```


