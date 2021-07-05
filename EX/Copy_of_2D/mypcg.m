function [x,iter,resvec] = mypcg(A,b,x0,maxit,tol,L)

    
    
    r = b - A*x0;
    p = L'\(L\r);
    rho = r.'*p;
    x = x0;
    err = tol + 1;
    tolB = tol*norm(b);
    resvec =[];
    resvec = [resvec,norm(r)];
    iter = 0;
    while (iter < maxit && err > tolB)
        
        
        z = A*p;
        alpha = rho/(z.'*p);
        x = x + alpha*p;
        r = r - alpha*z;
        g = (L')\(L\r);
        rho_new = r.'*g;
        beta = rho_new/rho;
        p = g + beta*p;
        
        %new step init
        rho = rho_new;
        err = norm(r);
        resvec = [resvec,err];
        iter = iter + 1;
            
            
    end
    
end

