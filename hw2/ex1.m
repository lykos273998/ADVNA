N = 201;
A = delsq(numgrid('S', N + 1));
n = size(A,1);
b = ones(n,1)/(n);
L = ichol(A);

tol = 1e-6;
maxit = 2000;

[x1, flag1, relres1, iter1, resvec1] = pcg(A,b,tol,maxit,L,L');
[x2, resvec2, it] = mypcg(A,b,tol,maxit,L);

semilogy(1:size(resvec1), resvec1,"g-o")
hold on 
semilogy(1:size(resvec2,2), resvec2, "b-*")
hold off
legend("MATLAB pcg - IC(0) ","mypcg - IC(0) ")


fprintf("|| sol_mypcg - sol_pcg || = %e \n", norm(x1-x2) )