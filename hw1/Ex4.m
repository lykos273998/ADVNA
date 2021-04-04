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

outfile = fopen("results_ex4.txt", "w");

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

u0_BWE_SOR = u0(x);
u0_CN_SOR = u0(x);



%intialization of the Laplacian matrix
L = L_1D(N);



for M = Time_subdivisions_BWE_CN
    ht = T/M;
    u_BWE_SOR = u0_BWE_SOR;
    u_CN_SOR = u0_CN_SOR;
    
    

    LHS_MAT_BWE = eye(N-1);
    LHS_MAT_BWE= LHS_MAT_BWE -  c * ht / (hx^2) * L;
    
    LHS_MAT_CN = eye(N-1);
    AUX_CN = 0.5 * c * ht / (hx^2) * L;
    LHS_MAT_CN= LHS_MAT_CN - AUX_CN;
    
    
    Tot_it_BWE_SOR = zeros(1,M);
    Tot_it_CN_SOR = zeros(1,M);
    
    %eigenvalues calculation and omega opt calculation
    zeta = c*ht/(hx^2);
    l1_cn = zeta/(1+zeta)*(1 - 2*(sin(pi/(2*N))^2));
    l2_cn = zeta/(1+zeta)*(1 - 2*(sin(pi*(N-1)/(2*N))^2));
    
    rhoJ_cn = max(abs(l1_cn), abs(l2_cn));
    omega_cn =  2/(1 + sqrt(1 - rhoJ_cn));
    
    l1_bwe = zeta/(1+2*zeta)*(2 - 4*(sin(pi/(2*N))^2));
    l2_bwe = zeta/(1+2*zeta)*(2 - 4*(sin(pi*(N-1)/(2*N))^2));
    
    rhoJ_bwe = max(abs(l1_bwe), abs(l2_bwe));
    omega_bwe =  2/(1 + sqrt(1 - rhoJ_bwe));
    
    
    
    for i = 1:M
        t = i*ht;
                
        tic();
        [u_BWE_SOR it_BWE_SOR] = BWE_SOR(u_BWE_SOR, x, t, LHS_MAT_BWE, c, hx, ht, omega_bwe); 
        elapsed_BWE_SOR = toc()*1000;
        
        tic();
        [u_CN_SOR it_CN_SOR] = CN_SOR(u_CN_SOR, x, t, LHS_MAT_CN, c, hx, ht, AUX_CN, omega_cn); 
        elapsed_CN_SOR = toc()*1000;
        
        Tot_it_BWE_SOR(i) = it_BWE_SOR;
        Tot_it_CN_SOR(i) = it_CN_SOR;
    end
   
    eBWE_SOR = norm(exact_solution - u_BWE_SOR);
    eCN_SOR = norm(exact_solution - u_CN_SOR);
    
    fprintf(outfile,"M = %d\n", M);
    fprintf(outfile,"--------------------------------------------------------------------------------------------\n");
    fprintf(outfile,"BWE SOR      Tot iter: %d\t Avg iter: %.2f\t Accuracy: %e\t CPU time(ms): %.2f\n", ...
        sum(Tot_it_BWE_SOR), mean(Tot_it_BWE_SOR), eBWE_SOR, elapsed_BWE_SOR);
    fprintf(outfile,"CN SOR       Tot iter: %d\t Avg iter: %.2f\t Accuracy: %e\t CPU time(ms): %.2f\n\n", ...
        sum(Tot_it_CN_SOR), mean(Tot_it_CN_SOR), eCN_SOR, elapsed_CN_SOR);
    
end
disp("You can find the results in `results_ex4.txt`");


function [u_m1, it] = CN_SOR(u_m, x, t, A, c, hx, ht, AUX_CN, omega)
    rhs = u_m + AUX_CN*u_m + 0.5*(f(x,t) + f(x, t + ht));
    tol = 1e-3 * ht;
    [u_m1, it, vdiff] = SOR(A, rhs, u_m, int16(1e5), tol,omega); 
end

function [u_m1, it] = BWE_SOR(u_m, x, t, A, c, hx, ht, omega)
    ff = ht * f(x, t + ht);
    rhs = u_m +  ff;
    tol = 1e-3 * ht;
    [u_m1, it, vdiff] = SOR(A, rhs, u_m, int16(1e5), tol, omega);
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