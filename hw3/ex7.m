Q = load('Q.dat');
Q = spconvert(Q);

Q = Q + Q' - sparse(diag(diag(Q)));

A = load('A.dat');
A = spconvert(A);

tol = 1e-12;
maxit = 500;


n = size(Q,1);
m = size(A,1);

H = @(x) [Q*x(1:n) + (A')*x(n+1:end); A*x(1:n)];

E = sparse(diag(diag(Q)));
E_1 = inv(E);

S = A*E_1*A';

p = amd(S);
S0 = S(p,p);
U = chol(S0);

M = @(x) [E\x(1:n); applyschur(x(n+1:end),p,U)];



x_exact = [1e-6*rand(n,1); rand(m,1)];
b = H(x_exact);
x0 = M(b);
[x,flag,relres,iter,resvec] = minres(H,b,tol,maxit,M,[],x0); %,

semilogy(resvec, '.-');
ylabel("||r_k||");
xlabel("iterations");
title("MINRES convergence profile");

function w = applyschur(v,p,U)
    n = size(U,1);
    z = v(p);
    u = U\(U'\z);
    pinv(p) = 1:n;
    w = u(pinv);
end

