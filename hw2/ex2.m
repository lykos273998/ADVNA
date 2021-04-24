Ns = [101,201,301,401];
%Ns = [4]
outfile = fopen('Ex2_results.txt','w');

fprintf(outfile, "                        |CG_________|______________PCG_________________\n");
fprintf(outfile, "                        |           |IC          IC(1e-2)    IC(1e-3)   \n");
fprintf(outfile, "N-1 	 n       k(A)   |Iter   CPU |Iter   CPU |Iter   CPU |Iter   CPU |\n");
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
    
    L = ichol(A);
    tic()
    [x2, flag2, relres2, iter2, resvec2] = pcg(A,b,tol,maxit,L,L');
    t2 = toc();
    iter_est = - log(err_red)/4*(sqrt(kA) + 1);
    its = [its,int32(iter_est)];
    
    %ICT tol 1e-2
    
    opts3.type = 'ict';
    opts3.droptol = 1e-2;
    L = ichol(A,opts3); 
    %whos L
    tic();
    [x3, flag3, relres3, iter3, resvec3] = pcg(A,b,tol,maxit, L, L');
    t3 = toc();
    
    %ICT tol 1e-3
    
    opts4.type = 'ict';
    opts4.droptol = 1e-3;
    L = ichol(A,opts4) ;
    tic();
    [x4, flag4, relres4, iter4, resvec4] = pcg(A,b,tol, maxit, L,L');
    t4 = toc();
    
    
    
                    %"N-1    n   |Iter      CPU|Iter      CPU|Iter      CPU|Iter      CPU|\n");
    fprintf(outfile, "%d \t%d \t %d\t|%d \t%.2f|%d \t%.2f|%d \t%.2f|%d \t%.2f|\n", ...
            N - 1, n, int32(kA), iter1, t1, iter2, t2,iter3, t3,iter4, t4);
end


disp("k estimated")
disp(its)
disp("! results in 'Ex2_results.txt' !")


