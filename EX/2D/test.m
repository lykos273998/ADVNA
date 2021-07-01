N=512;
A=delsq(numgrid('S',N+1));

n = size(A,1);

x0 = zeros(n,1);

x_exc = ones(n,1);
b = A*x_exc;

maxit = 10000;
tol = 1e-10;
tic();
[x, iter,resvec] = MG_2D_gen(A,b,x0,maxit,tol,N - 1,2);
t2 = toc()
tic();
[x2, iter2,resvec2] = SOR(A,b,x0,3,tol,1);
%xx = A\b;
t3 = toc()
iter
iter2

norm(x - x_exc)

semilogy(resvec)

x_p = linspace(0,1,15).';


rate = -log10(norm(x_exc - x)/norm(x_exc - x0))/iter









