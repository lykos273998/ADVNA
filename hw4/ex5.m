F = @(x) [x(1)^2 + x(2)^2 - 4;  x(1)*x(2) - 1];

Jac =@(x) [ 2*x(1) 2*x(2); x(2) x(1)];
x0 = [0; 1];
J0 = Jac(x0);
B0 = @(x) J0;



[x, iter, resvec] = quasi_newton(x0,F,B0,1e-8,100);

semilogy(resvec,'*-');
