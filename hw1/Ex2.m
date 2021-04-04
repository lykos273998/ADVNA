% Use the Jacobi and Gauss-Seidel methods to solve the linear system Ax = b,
% arising from the discretization of the Poisson equation in the unit square

%--- ex 2.1 ---
%Comparison between Gauss-Siedel and Jacobi
%GS is implemented as SOR with omega = 1


tol = 1e-8;
max_iter = 1000000;
disp("Exercise 2.1");
for N = [21,41,61,81]
    n_lowercase = (N - 1)*(N - 1);
    x0 = zeros(n_lowercase, 1);
    b = ones(n_lowercase, 1)/(N*N);
    A = delsq(numgrid('S', N + 1));
    
    [xj, IterJacobi, conv_profile_Jacobi] = jacobi(A,b,x0, max_iter, tol);
    [xGS, IterGS, conv_profile_GS] = SOR(A,b,x0, max_iter, tol, 1);
    
    figure
    semilogy(1:IterJacobi, conv_profile_Jacobi,'g')
    hold on
    semilogy(1:IterGS, conv_profile_GS,'b')
    hold off
    title(sprintf("N = %d, System dimension = %d",N, n_lowercase)) 
    legend('Jacobi','GS')
    
end


%--- ex 2.2 ---
%
disp("Exercise 2.2");
for N = [21,41,61,81]
    h = 1/N;
    rho_Hj = 2*(cos(h * pi / 2)^2) - 1;
    
    n_lowercase = (N - 1)*(N - 1);
    x0 = zeros(n_lowercase, 1);
    b = ones(n_lowercase, 1)/(N*N);
    A = delsq(numgrid('S', N + 1));
    
    [xj, IterJacobi, conv_profile_Jacobi] = jacobi(A,b,x0, max_iter, tol);
    
    x_true = A\b;
    true_error = norm(xj - x_true);
    estimated_error = conv_profile_Jacobi(IterJacobi) * rho_Hj / (1 - rho_Hj);
    disp(sprintf("System size = %d \t Rho_Hj = %f \t Estimated error %.7f \t True error %.7f",n_lowercase, rho_Hj, estimated_error, true_error));
end

%--- ex 2.3 ---
%Tridiagonal matrices arising from discretization of PDE are 2 cyclic 
%and consinstently oredered, so Youg-Varga theorem holds in this case
%so
disp("Exercise 2.3")
N = 41;
h = 1/N;
rho_Hj = 2*(cos(h * pi / 2)^2) - 1;
rho_GS = rho_Hj^2;

omega_opt = 2/(1 + sqrt(1 - rho_GS));
rho_SOR = omega_opt - 1;



n_lowercase = (N - 1)*(N - 1);
x0 = zeros(n_lowercase, 1);
b = ones(n_lowercase, 1)/(N*N);
A = delsq(numgrid('S', N + 1));

%using result obtained in exercise sheet 2
e0_est = norm(b/(h^2) - A*x0/(h^2))/18;
reduction = tol/e0_est;
p = - log10(reduction);
disp(sprintf("Reduction factor %f",p));
est_iter_J = int16(p / ( - log10(rho_Hj)));
est_iter_GS = int16(p / ( - log10(rho_GS)));
est_iter_SOR = int16(p / ( - log10(rho_SOR)));

[xj, IterJacobi, conv_profile_Jacobi] = jacobi(A,b,x0, max_iter, tol);
[xGS, IterGS, conv_profile_GS] = SOR(A,b,x0, max_iter, tol, 1);
[xSOR, IterSOR, conv_profile_SOR] = SOR(A,b,x0, max_iter, tol, omega_opt);

disp(sprintf("Jacobi Iterations -> estimated: %d \t done: %d", est_iter_J, IterJacobi));
disp(sprintf("GS Iterations     -> estimated: %d \t done: %d", est_iter_GS, IterGS));
disp(sprintf("SOR Iterations    -> estimated: %d \t done: %d", est_iter_SOR, IterSOR));


