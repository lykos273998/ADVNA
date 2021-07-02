function [x, iter,resvec] = MG_2D_gen(A,b,x0,maxit,tol,max_lvl)
    %2 level multigrid test
    iter = 0;
    exit_tol = tol*norm(b);
    x = x0;
    res = norm(b - A*x);
    resvec = [ res ];
    
    
    n = size(A,1); 
    
    
    for i = 2:max_lvl
        A_vec{i} = L1D(floor(n/(2^(i-1))) + 1);
    end
        
    nu1 = 5;
    nu2 = 5;
    
    nu_tot = 0;
        
    while iter < maxit && res > exit_tol
        iter = iter + 1;
        
        %first level
        [x, i, vdiff] = jacobi(A,b,x,nu1,tol);
        
        r_h = b - A*x;
        
        
        r_2h = f_to_c(r_h);
        
        
        %recursively call the V_cycle function
        
        e_2h = V_cycle_1D(r_2h, nu1, nu2, A_vec, tol, 2, max_lvl) ;    
      
        e_h = c_to_f(e_2h);
        
        x = x + e_h;
             
        [x, i, vdiff] = jacobi(A,b,x,nu2,tol);
        
        res = norm(b - A*x);
        resvec = [resvec,res];
        
    end
    
end

