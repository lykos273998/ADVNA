
F = @(x) [x(1)^2 + x(2)^2 - 4;  x(1)*x(2) - 1];

Jac =@(x) [ 2*x(1) 2*x(2); x(2) x(1)];

x0 = [0; 1];

[x, iter, resvec] = newton(x0,F,Jac,1e-8,10,1);

