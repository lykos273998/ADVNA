N=16;
A=delsq(numgrid('S',N+1));

B = delsq(numgrid('S', N/2 + 1));
size(B,1)
[P,R] = get_PR_fw(N/2 - 1); 

size(P);

B;

R*A*P;


full(R);

%gg = [1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 2 3 3 3 2 1 1 2 3 4 3 2 1 1 2 3 3 3 2 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 ].';

%gg2 = [2 2 2 2 4 2 2 2 2].';

%norm(P*gg2 - c_to_f_2d(gg2, N/2 - 1));



N=1024;
A=delsq(numgrid('S',N+1));

n = size(A,1);

x0 = zeros(n,1);

x_exc = ones(n,1);
b = A*x_exc;

maxit = 10000;
tol = 1e-10;
tic();
[x, iter,resvec] = MG_2D_gen(A,b,x0,maxit,tol,N,128,1);
t2 = toc()

iter



rate = -log10(norm(x_exc - x)/norm(x_exc - x0))/iter

norm(x_exc - x)















