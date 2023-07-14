println("--- Study periodic orbits: ---")

let
@unpack nodes, D1, ic, i1, i2 = p.params

Re_start = 350
p.params.Re = Re_start
h = 0.5

Res = [348, 346, 344, 342, 340, 338, 336, 334, 332, 330]
tperiod = 135 # approx period time
period = Int(135 / h)
timesteps = 4 * period 

Ψi =  1e-3 * randn((n + 1), (n + 1)) # slightly perturbed initial solution
sol, times = CF.timestepping(Ψi, p, h, timesteps; save = true, verbose = true)
filename = "$folderpo/orbit_$(Re_start)_start.jld" 
@save filename sol times

for Re in Res
    println("PO: Re = $Re")
    i = argmax(sol[period:end, ic, ic])
    Ψi = sol[i, :, :] # start with maximum psi value after one period

    p.params.Re = Re
    sol, times = CF.timestepping(Ψi, p, h, timesteps; save = true, verbose = true)

    filename = "$folderpo/orbit_$(Re).jld" 
    @save filename sol times
end
end
