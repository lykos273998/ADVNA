A = load('mat13041.rig');
A = spconvert(A);

n = size(A,1);

x_exact = zeros(n,1);

for i = 1:n
    x_exact(i) = 1/sqrt(i);
end

b = A*x_exact;

tol = 1e-10;
maxit = 550;
setup.type = 'crout';
setup.droptol = 0.1;
[L,U] = ilu(A,setup);
x0 = zeros(n,1);

[x1 ,iter1,resvec1,flag1] = myprecgmres(A,b,tol,maxit,x0,'L',L,U);

[x2 ,flag2,res2,iter2,resvec2] = gmres(A,b,maxit,tol,1,L,U);

[x3 ,iter3,resvec3,flag3] = myprecgmres(A,b,tol,maxit,x0,'R',L,U);

A_handle = @(x) A*(U\(L\x));
[x4 ,flag4,res4,iter4,resvec4] = gmres(A_handle,b,maxit,tol,1);
x4 = U\(L\x4);

[x5 ,iter5,resvec5,flag5] = myprecgmres(A,b,tol,maxit,x0,'S',L,U);

A_handle = @(x)L\(A*(U\x));
[x6 ,flag6,res6,iter6,resvec6] = gmres(A_handle,L\b,maxit,tol,1);
x6 = U\x6;

%disp("native")
%[x2,flag2,relres2,iter2,resvec2] = gmres(A,b,550,tol);
r1 = b - A*x1;
r2 = b - A*x2;
r3 = b - A*x3;
r4 = b - A*x4;
r5 = b - A*x5;
r6 = b - A*x6;

outfile = fopen('ex3_results.txt','w');
fprintf(outfile,'MyGMRES                                                  MATLAB GMRES\n');
fprintf(outfile,'Last res        True res         Transf res      iter 	 Last res       True res        Transf res     iter\n');
fprintf(outfile,'-------- LEFT Preconditioner\n');
fprintf(outfile,'%e \t %e \t %e \t %d \t %e \t %e \t %e \t %d \n', ...
    resvec1(end), norm(r1), norm(U\(L\r1)), iter1, resvec2(end),norm(r2), norm(U\(L\r2)), iter2(2));

fprintf(outfile,'-------- RIGHT Preconditioner\n');
fprintf(outfile,'%e \t %e \t %e \t %d \t %e \t %e \t %e \t %d \n', ...
    resvec3(end), norm(r3), norm(r3), iter3,  resvec4(end), norm(r4), norm(r4), iter4(2));
fprintf(outfile,'-------- SPLIT Preconditioner\n');
fprintf(outfile,'%e \t %e \t %e \t %d \t %e \t %e \t %e \t %d \n', ...
    resvec5(end), norm(r5), norm(L\r5), iter5, resvec6(end), norm(r6), norm(L\r6), iter6(2) );

%semilogy(resvec1,'.-')

