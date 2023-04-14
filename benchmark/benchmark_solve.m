NX = 64;
NY = 64;
Re = 100;

c = Fparam(NX,NY,Re);

psi = readmatrix('64x64_initial_guess_Re100.txt');
psii = psi(3:NX-1,3:NY-1);    

u = psii(:);
tolmax = 1e-10;
itmax = 100;

f = @() FbihCF(u,c);
timeit(f)

jac = @() FjacDF1(u,@FbihCF,c);
timeit(jac)

newton = @() Fnewtn(u,tolmax,itmax,@FbihCF,c);
timeit(newton)
