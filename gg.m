A = delsq(numgrid('S',41));
n = size(A,1);
x_exc = ones(n,1);

b = A*x_exc;

tol = 1e-10;

[x_k,iter,resvec,flag] = GMRESrestart(A,b,50,20,x0,tol)

semilogy(resvec)