A = load('mat13041.rig');
A = spconvert(A);

n = size(A,1);

x_exact = zeros(n,1);

for i = 1:n
    x_exact(i) = 1/sqrt(i);
end

b = A*x_exact;
%parameters
tol = 1e-12;
maxit = 600;
setup.type = 'crout';
setup.droptol = 0.01;
[L,U] = ilu(A,setup);
%to keep constant the total number of iterations
%maxit = restart*outer so...
outer = maxit/10;
tic();
[x1 ,flag1,res1,iter1,resvec1] = gmres(A,b,10,tol,outer,L,U);
t1 = toc();

outer = maxit/20;
tic();
[x2 ,flag2,res2,iter2,resvec2] = gmres(A,b,20,tol,outer,L,U);
t2 = toc();

outer = maxit/30;
tic();
[x3 ,flag3,res3,iter3,resvec3] = gmres(A,b,30,tol,outer,L,U);
t3 = toc();

outer = maxit/50;
tic();
[x4 ,flag4,res4,iter4,resvec4] = gmres(A,b,50,tol,outer,L,U);
t4 = toc();
outfile = fopen('ex4_results.txt','w');

fprintf(outfile,'             last res   \t CPU time \t Iterations\n');
fprintf(outfile,'restart: 10 %e \t %e \t %d \n', resvec1(end), t1, (iter1(1) - 1)*10 + iter1(2));
fprintf(outfile,'restart: 20 %e \t %e \t %d \n', resvec2(end), t2, (iter2(1) - 1)*20 + iter2(2));
fprintf(outfile,'restart: 30 %e \t %e \t %d \n', resvec3(end), t3, (iter3(1) - 1)*30 + iter3(2));
fprintf(outfile,'restart: 50 %e \t %e \t %d \n', resvec4(end), t4, (iter4(1) - 1)*50 + iter4(2));

semilogy(resvec1,'.-r')
hold on
semilogy(resvec2,'.-g')
hold on
semilogy(resvec3,'.-b')
hold on
semilogy(resvec4,'.-m')
hold on

legend('restart: 10','restart: 20','restart: 30','restart: 50');
title("GMRES `mat13041.rig` w. ILU (droptol 1e-2) ");
hold off
