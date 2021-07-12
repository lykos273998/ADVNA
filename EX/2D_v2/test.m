N=16;
A=delsq(numgrid('S',N+1));

B = delsq(numgrid('S', N/2 + 1));
size(B,1);
P = get_PR_fw(N/2 - 1); 
P2 = get_PR_fw_opt(N/2 - 1);

size(P);



gg = [1 1 1 1 1 1 1 1 2 2 2 2 2 1 1 2 3 3 3 2 1 1 2 3 4 3 2 1 1 2 3 3 3 2 1 1 2 2 2 2 2 1 1 1 1 1 1 1 1 ].';

gg2 = rand((N/2 - 1)^2,1) ;

norm(P*gg2 - P2*gg2)



N=1024;
A=delsq(numgrid('S',N+1));

n = size(A,1);

x0 = zeros(n,1);

x_exc = ones(n,1);
b = A*x_exc;

maxit = 10000;
tol = 1e-10;
tic();
[x, iter,resvec] = MG_2D_gen_v2(A,b,x0,maxit,tol,N,16,1);
t2 = toc()

iter



rate = -log10(norm(x_exc - x)/norm(x_exc - x0))/iter

norm(x_exc - x)















