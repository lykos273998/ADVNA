N = 49;
A = delsq(numgrid('S',N + 1));
n = size(A,1);
F = @(x) A * x - 0.1 * sin(x) - 5;
J = @(x) A - spdiags(0.1 * cos(x),0,n,n);


x0 = 1000 * sin((1:n)');

tol = 1e-12;
itmax = 100;
[xstar_dir, iter_dir, resvec_dir]  = newton(x0,F,J,tol,itmax,1);
[xstar_gmres, iter_gmres, resvec_gmres]  = newton(x0,F,J,tol,itmax,0,50,20,1e-8);

semilogy(resvec_dir, '.-')
hold on
semilogy(resvec_gmres, '.-')
legend('direct method', 'gmres')
ylabel("||F(x_k)||")
xlabel("iteration")
hold off