function [xstar, iter, resvec] = newton(x0,F,Jac, tol, itmax,lsol, gmres_restart, gmres_maxit, gmres_tol)
    exit_tol = tol*norm(F(x0));
    res = exit_tol + 1;
    resvec = [];
    iter = 0;
    xstar = x0;
    f = -F(x0);
    
    if nargin < 7
        gmres_restart = 50;
        gmres_tol = 1e-10;
        gmres_maxit = 50;
    end
    
    
    switch lsol
        case 1 
            while (res > exit_tol & iter < itmax)
                iter = iter + 1;
                
                s = Jac(xstar)\f;
                xstar = xstar + s;
                f = -F(xstar);                
                res = norm(f);
                
                resvec = [resvec,res];
            end
  
        case 0
            while res > exit_tol & iter < itmax
                iter = iter + 1;
                J = Jac(xstar);
                [L,U] = ilu(J);
                
                [s, it_gmres] = gmres(J, f, gmres_restart, gmres_tol, gmres_maxit, L, U, xstar);
                xstar = xstar + s;
                
                f = -F(xstar);                
                res = norm(f);
                
                resvec = [resvec,res];
            end
            
        otherwise
            disp("not implemented case, use 1: direct solver, 0: iterative solver (GMRES)")
    end
    
    
    
end

