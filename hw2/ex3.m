Ns = [101,201,301,401];
%Ns = [4]
outfile = fopen('Ex3_results.txt','w');

fprintf(outfile, "                        |CG_________|PCG________|\n");
fprintf(outfile, "                        |           |Jacobi\n");
fprintf(outfile, "N-1 	 n       k(A)   |Iter   CPU |Iter   CPU |\n");
ks = [];
its = [];

for N = Ns
    A = delsq(numgrid('S',N+1));
    n = size(A,1);
    disp(sprintf("running: n = %d",n))
    tol = 1e-8;
    maxit = 2000;
    
    x_exact = zeros(n,1);
    for i = 1:n
        x_exact(i,1) = 1/sqrt(i);
    end
    
    b = A*x_exact;
    
    tic();
    %no preconditioning
    [x1, flag1, relres1, iter1, resvec1] = pcg(A,b,tol,maxit);
    t1 = toc();
    % IC(0)
    x_diff = x1 - x_exact;
    
   
    
    
    %kA = ( 2 * iter1 / log(ratio/2) + 1) ^2;
    kA = (4/(pi^2) * (N^2));
    ks = [ks, kA];
    err_red = sqrt(x_diff.'*A*x_diff)/sqrt(x_exact.'*A*x_exact);
    
    L = sparse(diag(diag(A)));
    tic()
    [x2, flag2, relres2, iter2, resvec2] = pcg(A,b,tol,maxit,L);
    t2 = toc();
    iter_est = - log(err_red)/4*(sqrt(kA) + 1);
    its = [its,int32(iter_est)];
    
   
    
    
    
                    %"N-1    n   |Iter      CPU|Iter      CPU|\n");
    fprintf(outfile, "%d \t%d \t %d\t|%d \t%.2f|%d \t%.2f|\n", ...
            N - 1, n, int32(kA), iter1, t1, iter2, t2);
end


disp("k estimated")
disp(its)
disp("! results in 'Ex3_results.txt' !")
fclose('all');


