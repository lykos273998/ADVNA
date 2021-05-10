h = 6.25e-3;
N = (1 + 1/h);
A = delsq(numgrid('S', N + 1));

n = size(A,1);

b = ones(n,1)*(h^2);
tol_stationary = 1e-10;
tol_cg = 1e-8;

x0 = zeros(n,1);
maxiter = n;

x_exact = A\b;
x_ex_norm = norm(x_exact);

outfile = fopen("Ex5_results.txt","w");
fprintf(outfile, "Method  \t Iterations \tCPU \t rel err \t optional info\n");

%a jacobi
disp("Running Jacobi")
tic();
[x1, iter1, vdiff1] = jacobi(A,b,x0,maxiter,tol_stationary);
t = toc();
fprintf(outfile,"Jacobi   \t: %d \t %.2f \t %e\n", iter1, t, norm(x1 - x_exact)/x_ex_norm);


%b Gauss seidel
disp("Running GS")
tic();
[x2, iter2, vdiff2] = SOR(A,b,x0,maxiter,tol_stationary,1);
t = toc();
fprintf(outfile,"Gauss-Seidel \t: %d \t %.2f \t %e\n", iter2, t, norm(x2 - x_exact)/x_ex_norm);

%c SOR
disp("Running SOR opt")
rho_j = max(abs(1 - 2 * (sin(pi/2 * h))^2),abs(1 - 2 * (sin(pi/2 * (1-h)))^2));

omega_opt = 2/(1 + sqrt(1 - rho_j^2));


tic();
[x3, iter3, vdiff3] = SOR(A,b,x0,maxiter,tol_stationary,omega_opt);
t = toc();
fprintf(outfile,"SOR (opt) \t: %d \t %.2f \t %e \t omega opt: %.2f rho_j %.4f rho_omega %.4f\n", iter3, t, norm(x3 - x_exact)/x_ex_norm, omega_opt, rho_j, omega_opt - 1);

%d SOR

disp("Running SOR * 1.01")
omega = omega_opt*1.01;

tic();
[x4, iter4, vdiff4] = SOR(A,b,x0,maxiter,tol_stationary,omega);
t = toc();

fprintf(outfile,"SOR (2) \t: %d \t %.2f \t %e \t omega: %.2f\n", iter4, t, norm(x4 - x_exact)/x_ex_norm, omega);

%e SOR


disp("Running SOR / 1.01")
omega = omega_opt/1.01;

tic();
[x5, iter5, vdiff5] = SOR(A,b,x0,maxiter,tol_stationary,omega);
t = toc();
fprintf(outfile,"SOR (3) \t: %d \t %.2f \t %e \t omega: %.2f\n", iter5, t, norm(x5 - x_exact)/x_ex_norm, omega);

%f PCG

disp("Running CG")

tic();
[x6, flag6, relres6, iter6, resvec6] = pcg(A,b,tol_cg,maxiter);
t = toc();
fprintf(outfile,"Conj. Grad. : %d \t %.2f \t %e \n", iter6, t, norm(x6 - x_exact)/x_ex_norm);

disp("Running PCG IC(0)")
%g PCG IC(0)
L = ichol(A);


tic();
[x7, flag7, relres7, iter7, resvec7] = pcg(A,b,tol_cg,maxiter, L, L');
t = toc();

fprintf(outfile,"PCG IC(0)  \t: %d \t %.2f \t %e \n", iter7, t, norm(x7 - x_exact)/x_ex_norm);
%h PCG ICT 1e-2
disp("Running PCG ICT")
opts.type = 'ict';
opts.droptol = 1e-2;
L = ichol(A, opts);
tic();
[x8, flag8, relres8, iter8, resvec8] = pcg(A,b,tol_cg,maxiter, L, L');
t = toc();

fprintf(outfile,"PCG ICT \t: %d \t %.2f \t %e \t tolerance: %e\n", iter8, t, norm(x7 - x_exact)/x_ex_norm, opts.droptol);

figure()
semilogy(1:iter3, vdiff3,'.-', 1:iter4, vdiff4,'.-', 1:iter5, vdiff5,'.-')
legend("SOR (\omega_{opt})", "SOR (\omega_{opt} * 1.01)", "SOR (\omega_{opt} / 1.01)")
xlabel("Iterations")
ylabel("||x_{k+1} - x_k||")
figure()
semilogy(1:(iter6 + 1),resvec6,'.-',1:(iter7 + 1),resvec7,'.-',1:(iter8 + 1),resvec8,'.-')
legend("CG", "PCG IC(0)", "PCG ICT (tol 1e-2)")
xlabel("Iterations")
ylabel("||r_k||")
fclose('all');










function [x, iter, vdiff] = SOR(A,b,x0,max_iter,tol,omega)
    %SOR solver function
    
    %decomposition of A
    L = tril(A,-1);
    U = triu(A,1);
    D = diag(diag(A));
    M = omega*L + D;
    N = (1 - omega)*D - omega*U;
    omega_b = omega*b;
    
    %initialization of auxiliary vectors
    x_old = x0;
    x_new = x0;
    err = tol + 1;
    vdiff = [];
    
    for iter = 1:max_iter
        
        
        x_temp = N*x_old + omega_b;
        x_new = M\x_temp;
        
        err = norm(x_new - x_old);
        x_old = x_new;
        vdiff = [vdiff, err];
        if err < tol
            break
        end
        
    end
    
    x = x_new;
    
end


function [x,iter,vdiff] = jacobi(A, b, x0, max_iter, tolerance)
    L = tril(A,-1);
    U = triu(A,1);
    D = diag(diag(A));
    M = D;
    N = -(L+U);
    err = tolerance + 1;
    x_old = x0;
    x_new = x0;
    vdiff = [];
    for iter = 1:max_iter
        
        x_temp = N * x_old + b;
        x_new = M\x_temp;
        
        err = norm(x_new - x_old);
        x_old = x_new;
        vdiff = [vdiff, err];
        if err < tolerance
            break
        end
        
    end
    
    x = x_new;
    
end


