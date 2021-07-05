function [x, iter,resvec] = MG_2D_gen(A,b,x0,maxit,tol,nx,min_finest_grid, verbose)
    %general multi grid V-cycle scheme
    iter = 0;
    exit_tol = tol*norm(b);
    x = x0;
    res = norm(b - A*x);
    resvec = [ res ];
    
    %by default solve on 2 grids
    max_lvl = 2;
    
    nx_aux = (nx + 1)/4;
    
    while nx_aux > min_finest_grid
        nx_aux = nx_aux/2;
        max_lvl = max_lvl + 1;
    end
    
    n = size(A,1); 
    
    nx_aux = nx - 1;
    
    A_vec{1} = A;
    %calculation of the all discrete laplacians needed 
    for i = 2:max_lvl
        nx_aux = floor(nx_aux/2);
        [P_vec{i}, R_vec{i}] = get_PR_fw(nx_aux);
        A_vec{i} = R_vec{i}*A_vec{i - 1}*P_vec{i};
    end
    
  
    nu1 = 5;
    nu2 = 3;
    
    if verbose
        fprintf("system size %d \n", size(A,1));
        fprintf("mesh subdivisions %d \n", nx);
        fprintf("solving exactly at level %d with mesh subd %d \n", max_lvl, nx/(2^(max_lvl - 1)));
        
    end
    
    while iter < maxit && res > exit_tol
        iter = iter + 1;
        
        %first level
        [x, i, vdiff] = SOR(A,b,x,nu1,tol,1);
        
        r_h = b - A*x;
      
        r_2h = R_vec{2} * r_h;
        
        %recursively call the V_cycle function
        
        e_2h = V_cycle(r_2h, nu1, nu2, A_vec,P_vec, R_vec, tol, 2, max_lvl) ;    
      
        e_h = P_vec{2} * e_2h;
        
        x = x + e_h;
        
      
        [x, i, vdiff] = SOR(A,b,x,nu2,tol,1);
        
        res = norm(b - A*x);
        resvec = [resvec,res];
        
        
    end
    
end

