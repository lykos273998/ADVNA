function [xstar, iter, resvec] = inexact_newton(x0,F,Jac, tol, itmax, eta_max, lambda,gmres_restart, gmres_maxit)
    exit_tol = tol*norm(F(x0));
    res = exit_tol + 1;
    resvec = [];
    iter = 0;
    xstar = x0;
    f = -F(x0);
    
    if nargin < 8
        gmres_restart = 50;
        gmres_maxit = 50;
    end
    
    eta_k = @(f_k, f_k_1) min(eta_max, min(lambda*(f_k^2)/(f_k_1^2));
    
    
    %first iteration
    iter = iter + 1;
    J = Jac(xstar);
    [L,U] = ilu(J);

    [s, it_gmres] = gmres(J, f, gmres_restart, eta_max, gmres_maxit, L, U, xstar);
    xstar = xstar + s;

    f = -F(xstar);                
    res = norm(f);

    resvec = [resvec,res];
    


    while res > exit_tol & iter < itmax
        iter = iter + 1;
        J = Jac(xstar);
        [L,U] = ilu(J);
        eta = eta_k(resvec(end),resvec(end-1));
        [s, it_gmres] = gmres(J, f, gmres_restart, eta, gmres_maxit, L, U, xstar);
        xstar = xstar + s;

        f = -F(xstar);                
        res = norm(f);

        resvec = [resvec,res];
    end

    
end

