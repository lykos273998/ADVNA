A = load('ML_laplace.mtx');
A = spconvert(A);
n = size(A,1);
restart = 50;

d_tols = [2e-2, 1e-2, 3e-3, 1e-3, 1e-4, 1e-5];
outfile = fopen('ex5_results.txt', 'w');
fprintf(outfile,'drop tolerance 	 iter    tprec  	 tsol        tTOT        density\n');
tol = 1e-10;
x_exact = ones(n,1);
b = A*x_exact;

for droptol = d_tols
    setup.type = 'crout';
    setup.droptol = droptol;
    
    tic();
    [L,U] = ilu(A,setup);
    tprec = toc();
    
    tic()
    [x ,flag,res,iter,resvec] = gmres(A,b,restart,tol,50,L,U);
    tsol = toc();
    rho = (nnz(L) + nnz(U) - n)/nnz(A);
    fprintf(outfile,'%e \t %d \t %.3fs \t %.3fs \t %.3fs \t %.3f \n',...
        droptol, (iter(1)-1)*restart + iter(2), tprec, tsol, tprec+tsol, rho);
    semilogy(resvec,'.-')
    hold on 
end
legend("2e-2", "1e-2", "3e-3", "1e-3", "1e-4", "1e-5")
ylabel("||r_k||")
xlabel("iteration")
title("GMRES w. ILU different drop tolerances")
