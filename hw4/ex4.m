N = 50;
A = delsq(numgrid('S',N + 1));
n = size(A,1);

[u0,lambda0] = eigs(A,1,'sm');
u0 = u0 + ones(n,1)*1e-2;
lambda0 = lambda0 + 1e-2;

F = @(x) [(A - speye(n)*x(end))*x(1:end-1); x(1:end-1)'*A*x(1:end-1) - x(end)];
J = @(x) [ A - speye(n)*x(end) A*x(1:end-1); (0.5*A*x(1:end-1))' -1];
itmax = 20;
tol = 1e-12;
[xstar_gmres, iter_gmres, resvec_gmres]  = newton([u0;lambda0],F,J,tol,itmax,0,50,20,1e-8);
semilogy(resvec_gmres,'.-')
[u0,lambda0] = eigs(A,1,'sm');
%xstar_gmres(end) - lambda0

uu = J([u0;lambda0]);
hh= F([u0;lambda0]);