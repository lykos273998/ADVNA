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
