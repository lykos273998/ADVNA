Ns = [101,201,301,401];
%Ns = [4]
outfile = fopen('Ex4_results.txt','w');

fprintf(outfile, "            |CG_________|__________PCG__________\n");
fprintf(outfile, "            |           |Jacobi     IC(0)       \n");
fprintf(outfile, "N-1 n       |Iter   CPU |Iter   CPU |Iter   CPU |  \n");
ks = [];
its = [];

for N = [100]
    A = gallery('wathen', 100, 100);
    n = size(A,1);
    disp(sprintf("running: n = %d",n))
    tol = 1e-8;
    maxit = 2000;
    
    x_exact = rand(n,1);   
    b = A*x_exact;
    
    tic();
    %no preconditioning
    [x1, flag1, relres1, iter1, resvec1] = pcg(A,b,tol,maxit);
    t1 = toc();
    % IC(0)
    
    
    %jacobi
    L = sparse(diag(diag(A)));
    tic()
    [x2, flag2, relres2, iter2, resvec2] = pcg(A,b,tol,maxit,L,L');
    t2 = toc();
    iter_est = - log(err_red)/4*(sqrt(kA) + 1);
    its = [its,int32(iter_est)];
    
    %IC0
    
    L = ichol(A); 
    %whos L
    tic();
    [x3, flag3, relres3, iter3, resvec3] = pcg(A,b,tol,maxit, L, L');
    t3 = toc();
    
    %ICT tol 1e-3
    
    
    
                    %"N-1    n   |Iter      CPU|Iter      CPU|Iter      CPU|Iter      CPU|\n");
    fprintf(outfile, "%d \t%d \t|%d \t%.2f|%d \t%.2f|%d \t%.2f|\n", ...
            N - 1, n, iter1, t1, iter2, t2,iter3, t3);
end

disp("! results in 'Ex4_results.txt' !")


