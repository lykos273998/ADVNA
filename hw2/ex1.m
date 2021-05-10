N = 101;
A = delsq(numgrid('S', N + 1));
n = size(A,1);
b = ones(n,1)/(n);
L = ichol(A);

tol = 1e-6;
maxit = 2000;

[x1, flag1, relres1, iter1, resvec1] = pcg(A,b,tol,maxit,L,L');
[x2, iter2, resvec2] = mypcg(A,b,tol,maxit,L);

semilogy(1:(iter1 + 1), resvec1,"g-o")
hold on 
semilogy(1:(iter2 + 1), resvec2, "b-*")
hold off
legend("MATLAB pcg - IC(0) ","mypcg - IC(0) ")
title(sprintf("N = %d, System size = %d, tolerance = %.2e", N, n, tol));
xlabel('Iterations');
ylabel('||r_k||');

fprintf("|| sol_mypcg - sol_pcg || = %4.2e \n", norm(x1-x2) )
fprintf("|| resvec_pcg - resvec_mypcg || = %4.2e \n", norm(x1-x2) )
fprintf("pcg iteration %d resvec last %e\n", iter1, resvec1(iter1));
fprintf("MYpcg iteration %d resvec last %e\n", iter2, resvec2(iter2));

