function [xstar, iter, resvec] = newton_bt(x0,F,Jac, tol, itmax,lsol, tau,  gmres_restart, gmres_tol, gmres_maxit)
    exit_tol = tol*norm(F(x0));
    res = exit_tol + 1;
    resvec = [];
    iter = 0;
    xstar = x0;
    f = -F(x0);
    resvec = [resvec, norm(f)]
    if nargin < 8
        gmres_restart = 50;
        gmres_tol = 1e-10;
        gmres_maxit = 50;
    end
    
    
    switch lsol
        case 1 
            while (res > exit_tol & iter < itmax)
                iter = iter + 1;
                
                s = Jac(xstar)\f;
                alpha = 1;
                
                while norm(f) >= resvec(end)
                    x_candidate = xstar + alpha*s;
                    f = -F(x_candidate);
                    alpha = alpha*tau;
                end
                
                xstar = x_candidate
                                
                res = norm(f);
                
                resvec = [resvec,res];
            end
  
        case 0
            while res > exit_tol & iter < itmax
                iter = iter + 1;
                
                s = gmres(Jac(xstar), f, gmres_restart, gmres_tol, gmres_maxit);
                alpha = 1;
                
                while norm(f) >= resvec(end)
                    x_candidate = xstar + alpha*s;
                    f = -F(x_candidate);
                    alpha = alpha*tau;
                end
                
                xstar = x_candidate              
                res = norm(f);
                
                resvec = [resvec,res];
            end
            
        otherwise
            disp("not implemented case, use 1: direct solver, 0: iterative solver (GMRES)")
    end
    
    
    
end

