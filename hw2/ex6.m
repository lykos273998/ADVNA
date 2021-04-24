%matrix loading
B = load('apache1.mat');
matrix = cell2mat(struct2cell(B));
A = matrix.A;
x_exact = ones(size(A,1),1);
b = A*x_exact;
tol = 1e-8;
maxiter = 5000;
[x1, flag1, relres1, iter1, resvec1] = pcg(A,b,tol,maxiter);

L = ichol(A);

[x2, flag2, relres2, iter2, resvec2] = pcg(A,b,tol,maxiter, L, L');

b_norm = norm(b);
semilogy(1:(iter1 + 1), resvec1/b_norm ,1:(iter2 + 1), resvec2/b_norm);
legend("CG","PCG IC(0)");
xlabel("iterations")
ylabel("||r_k||/||b||")