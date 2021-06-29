N=1024;
%A=delsq(linspace(0,1,N+1));

A = L1D(N);
n = size(A,1);

x0 = zeros(n,1);

x_exc = ones(n,1);
b = A*x_exc;

maxit = 10000;
tol = 1e-10;
tic();
[x, iter,resvec] = MG_2L_1D(A,b,x0,maxit,tol);
t2 = toc()
tic();
[x2, iter2,resvec2] = SOR(A,b,x0,maxit,tol,1);
%xx = A\b;
t3 = toc()

iter
iter2
norm(x - x_exc)

%semilogy(resvec)

x_p = linspace(0,1,15).';
%x_pc = linspace(0,1,7).';
%x_pf = linspace(0,1,31).';


