function [x, iter,resvec] = MG_2D_gen(A,b,x0,maxit,tol,nx,max_lvl)
    %general multi grid V-cycle scheme
    iter = 0;
    exit_tol = tol*norm(b);
    x = x0;
    res = norm(b - A*x);
    resvec = [ res ];
    
    
    n = size(A,1); 
    
    nx_2 = nx;
    %calculation of the all discrete laplacians needed 
    for i = 2:max_lvl
        nx_2 = floor(nx_2/2);
        A_vec{i} = delsq(numgrid('S',nx_2 + 2));
    end
    
    
    
    
   
    
    nu1 = 5;
    nu2 = 3;
    
    nu_tot = 0;
    
    while iter < maxit && res > exit_tol
        iter = iter + 1;
        
        %first level
        [x, i, vdiff] = SOR(A,b,x,nu1,tol,1);
        
        r_h = b - A*x;
      
        r_2h = f_to_c_2d(r_h,floor(nx/2));
        
        %recursively call the V_cycle function
        
        e_2h = V_cycle(r_2h, floor(nx/2), nu1, nu2, A_vec, tol, 2, max_lvl) ;    
      
        e_h = c_to_f_2d(e_2h,floor(nx/2));
        
        x = x + e_h;
        
      
        [x, i, vdiff] = SOR(A,b,x,nu2,tol,1);
        
        res = norm(b - A*x);
        resvec = [resvec,res];
        
        
    end
    
end

