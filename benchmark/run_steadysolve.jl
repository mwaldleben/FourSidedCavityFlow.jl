n = 64
dim = n+1
nodes, D1, D2, D4 = diffchebyshev(n)

Re = 100

k0 = 10
bcfunc(x) = ((exp(k0*(x-1)) - 1) * (exp(-k0*(x+1)) - 1))^2

bcleft = bcfunc.(nodes)
bcright = -bcfunc.(nodes)
bctop = bcfunc.(nodes)
bcbottom = -bcfunc.(nodes)

Minv = constructBCmatrix(D1) 
h1 = similar(bctop[3:n-1])
h2 = similar(bctop[3:n-1])
k1 = similar(bctop[2:n])
k2 = similar(bctop[2:n])

Ψ = zeros(n+1, n+1)
fΨ = similar(Ψ)
laplΨ = similar(Ψ)
biharmΨ = similar(Ψ)
DΨ = similar(Ψ)
ΨD = similar(Ψ)
lDΨ = similar(Ψ)
lΨD = similar(Ψ)
nonlinΨ = similar(Ψ)

Ψ0i = 1e-3*randn(n-3,n-3)
u0 = vec(Ψ0i)
fu = similar(u0) 

p = CavityCache{Float64}(fΨ, Ψ, n, Re, D1, D2, D4, bcleft, bcright, bctop, bcbottom, Minv, h1, h2, k1, k2, DΨ, ΨD, lDΨ, lΨD, biharmΨ, laplΨ, nonlinΨ)

# Test custom solve
fnewton(u) = f(u,p) 
fnewton!(fu,u) = f!(fu,u,p) 
fnewton!(fu, u0)
sol, iter, tol = @time newton_custom!(fnewton!, u0)

# Plot
# Ψsol = zeros(n+1, n+1)
# Ψsol[3:n-1,3:n-1][:] = sol
# constructΨ!(Ψsol, n, D1, Minv, bcleft, bcright, bctop, bcbottom, h1, h2, k1, k2)
# contourf(reverse(nodes), reverse(nodes), Ψ')
