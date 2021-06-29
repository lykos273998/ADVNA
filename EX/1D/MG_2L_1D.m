function [x, iter,resvec] = MG_2L_1D(A,b,x0,maxit,tol)
    %2 level multigrid test
    iter = 0;
    exit_tol = tol*norm(b);
    x = x0;
    res = norm(b - A*x);
    resvec = [ res ];
    
    n = size(A,1);
    A2 = L1D(floor(n/2) + 1);
    A4 = L1D(floor(n/4) + 1);
   
    
    nu1 = 5;
    nu2 = 3;
    
    r_2h = zeros(floor(n/2),1);
    e_2h_0 = zeros(floor(n/2),1);
    e_4h_0 = zeros(floor(n/4),1);
    nu_tot = 0;
    
    while iter < maxit && res > exit_tol
        iter = iter + 1;
        
        %first level
        [x, i, vdiff] = SOR(A,b,x,nu1,tol,1);
        
        r_h = b - A*x;
      
        r_2h = f_to_c(r_h);
        
        %second level        
        [e_2h, i, vdiff] = SOR(A2,r_2h,e_2h_0,nu1,tol,1);
        
        
        r_4h = f_to_c(r_2h - A2*e_2h);
        
        %third level
        %[e_4h, i, vdiff] = SOR(A4,r_4h,e_4h_0,nu2,tol,1);
        e_4h = A4\r_4h;
        
        e_2h = e_2h + c_to_f(e_4h);
        
        [e_2h, i, vdiff] = SOR(A2,r_2h,e_2h,nu2,tol,1);
      
        e_h = c_to_f(e_2h);
        
        x = x + e_h;
      
        [x, i, vdiff] = SOR(A,b,x,nu2,tol,1);
        
        res = norm(b - A*x);
        resvec = [resvec,res];
        
        
    end
    
end

