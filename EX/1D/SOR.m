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

