N = 50;
A = delsq(numgrid('S',N + 1));
n = size(A,1);

[u0,lambda0] = eigs(A,1,'sm');

u0 = u0 + ones(n,1)*1e-2;
lambda0 = lambda0 + 1e-2;

F = @(x) [A*x(1:end-1) - x(end)*x(1:end-1); x(1:end-1).'*x(1:end-1) - 1];
J = @(x) [ A - speye(n)*x(end) -x(1:end-1); (2*x(1:end-1))' 0];
itmax = 20;
tol = 1e-12;
[xstar, iter, resvec]  = newton([u0;lambda0],F,J,tol,itmax,1);

semilogy(resvec,'.-')
ylabel("||F(x_k)||")
xlabel("iteration")

[uu0,ll0] = eigs(A,1,'sm');

disp(norm(xstar(1:end-1) - uu0))
disp(abs(ll0 - xstar(end)))
