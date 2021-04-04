%exercise 3

%/!\ NOTE /!\
%The code is not optimized for performance, the aim is to make it clear and
%understandable, there are lots of quantities computed multiple times that
%can be computed once and then reused
%/!\ NOTE /!\
%in the case considered in the exercise the values of u at the boundaries
%of the spatial domain are fixed to be 0

%parameters
T = 1/30;
c = 1;
N = 160;

outfile = fopen("results_ex3.txt", "w");

%setting up the time subdivions to test
Time_subdivisions_FWE = [600, 1200, 2400, 4800];
Time_subdivisions_BWE_CN = [12, 60, 120, 600];
exact_solution = load('truesol.txt');
exact_solution = exact_solution(2:N);


%slicing the domain cutting of the boundary
x = linspace(0,1, N + 1);
x = x(2:N);
x = x.';

hx = 1/N;

%initialization of initial states
u0_FWE = u0(x);
u0_BWE_J = u0(x);
u0_BWE_GS = u0(x);
u0_CN_GS = u0(x);



%intialization of the Laplacian matrix
L = L_1D(N);
ht = 0.1;
%Forward Euler testing:

for M = Time_subdivisions_FWE
    ht = T/M;
    u_FWE = u0_FWE;
    for i = 1:M
        t = i*ht;
        u_FWE = FWE(u_FWE, x, t, L, c, hx, ht);   
    end
    norm(exact_solution - u_FWE);
end


for M = Time_subdivisions_BWE_CN
    ht = T/M;
    u_BWE_J = u0_BWE_J;
    u_BWE_GS = u0_BWE_GS;
    u_CN_GS = u0_CN_GS;
    
    

    LHS_MAT_BWE = eye(N-1);
    LHS_MAT_BWE= LHS_MAT_BWE -  c * ht / (hx^2) * L;
    
    LHS_MAT_CN = eye(N-1);
    AUX_CN = 0.5 * c * ht / (hx^2) * L;
    LHS_MAT_CN= LHS_MAT_CN - AUX_CN;
    
    Tot_it_BWE_J = zeros(1,M);
    Tot_it_BWE_GS = zeros(1,M);
    Tot_it_CN_GS = zeros(1,M);
    
    for i = 1:M
        t = i*ht;
        
        tic();
        [u_BWE_J, it_BWE_J] = BWE_Jacobi(u_BWE_J, x, t, LHS_MAT_BWE, c, hx, ht);
        elapsed_BWE_J = toc()*1000;
        Tot_it_BWE_J(i) = it_BWE_J;
        
    end
     
    for i = 1:M
        t = i*ht;
        
        tic();
        [u_BWE_GS it_BWE_GS] = BWE_GS(u_BWE_GS, x, t, LHS_MAT_BWE, c, hx, ht); 
        elapsed_BWE_GS = toc()*1000;
        Tot_it_BWE_GS(i) = it_BWE_GS;
      
    end
        
    for i = 1:M
        t = i*ht;    
        
        tic();
        [u_CN_GS it_CN_GS] = CN_GS(u_CN_GS, x, t, LHS_MAT_CN, c, hx, ht, AUX_CN); 
        elapsed_CN_GS = toc()*1000;
        Tot_it_CN_GS(i) = it_CN_GS;
        
    end
       
    eBWE_J = norm(exact_solution - u_BWE_J);
    eBWE_GS = norm(exact_solution - u_BWE_GS);
    eCN_GS = norm(exact_solution - u_CN_GS);
    
    fprintf(outfile,"M = %d\n", M);
    fprintf(outfile,"--------------------------------------------------------------------------------------------\n");
    fprintf(outfile,"BWE Jacobi  Tot iter: %d\t Avg iter: %.2f\t Accuracy: %e\t CPU time(ms): %.2f\n", ...
        sum(Tot_it_BWE_J), mean(Tot_it_BWE_J), eBWE_J, elapsed_BWE_J);
    fprintf(outfile,"BWE GS      Tot iter: %d\t Avg iter: %.2f\t Accuracy: %e\t CPU time(ms): %.2f\n", ...
        sum(Tot_it_BWE_GS), mean(Tot_it_BWE_GS), eBWE_GS, elapsed_BWE_GS);
    fprintf(outfile,"CN GS       Tot iter: %d\t Avg iter: %.2f\t Accuracy: %e\t CPU time(ms): %.2f\n\n", ...
        sum(Tot_it_CN_GS), mean(Tot_it_CN_GS), eCN_GS, elapsed_CN_GS);
    
end
disp("You can find the results in `results_ex3.txt`");

function [u_m1, it] = CN_GS(u_m, x, t, A, c, hx, ht, AUX_CN)
    rhs = u_m + AUX_CN*u_m + 0.5*(f(x,t) + f(x, t + ht));
    tol = 1e-3 * ht;
    [u_m1, it, vdiff] = SOR(A, rhs, u_m, int16(1e5), tol,1); 
end

function [u_m1, it] = BWE_GS(u_m, x, t, A, c, hx, ht)
    ff = ht * f(x, t + ht);
    rhs = u_m +  ff;
    tol = 1e-3 * ht;
    [u_m1, it, vdiff] = SOR(A, rhs, u_m, int16(1e5), tol, 1);
end

function [u_m1, it] = BWE_Jacobi(u_m, x, t, A, c, hx, ht)
    ff = ht * f(x, t + ht);
    rhs = u_m +  ff;
    tol = 1e-3 * ht;
    [u_m1, it, vdiff] = jacobi(A, rhs, u_m, int16(1e5), tol);
end

function u_m1 = FWE(u_m, x, t, L, c, hx, ht)
    u_m1 = u_m + c * ht / (hx^2) * L * u_m + ht * f(x,t);
end

function L = L_1D(N)
    u = - 2 * ones(N-1,1);
    w = ones(N-2,1);
    L = sparse(diag(u) + diag(w,-1) + diag(w,1));
end



function u = u0(x,t)
    u = (x - x.*x) .* (x.*x + sin(2*pi*x));
end

function ff = f(x,t)
    ff = zeros(length(x),1);
end 