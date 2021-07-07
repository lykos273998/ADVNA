function [x, iter,resvec,t1,t2] = MG_2D_gen(A,b,x0,maxit,tol, min_mesh, verbose )
    %general multi grid V-cycle scheme
    iter = 0;
    exit_tol = tol*norm(b);
    x = x0;
    res = norm(b - A*x);
    resvec = [ res ];
    
    
    n = size(A,1); 
    
    nx_aux = (n + 1)/2;
    
    A_vec{1} = A;
    
    max_lvl = 2;
    
    tic();
    while nx_aux > min_mesh
        nx_aux = nx_aux/2;
        max_lvl = max_lvl + 1;
    end
        
    
    nx_aux = n;
    for i = 2:max_lvl
        nx_aux = floor(nx_aux/2);
        P_vec{i} = get_PR_opt(nx_aux);
        A_vec{i} =( P_vec{i}.' * A_vec{i-1} * P_vec{i} )/(2*4);
    end
    
    t1 = toc();
        
    nu1 = 5;
    nu2 = 3;
    
    nu_tot = 0;
       
    tic();
    while iter < maxit && res > exit_tol
        iter = iter + 1;
        
        %first level
        [x, i, vdiff] = jacobi(A,b,x,nu1,tol);
        
        r_h = b - A*x;
        
        
        r_2h = (P_vec{2}.' * r_h)/4;
        
        
        %recursively call the V_cycle function
        
        e_2h = V_cycle_1D(r_2h, nu1, nu2, A_vec,P_vec, tol, 2, max_lvl) ;    
      
        e_h = (P_vec{2} * e_2h)/2;
        
        x = x + e_h;
             
        [x, i, vdiff] = jacobi(A,b,x,nu2,tol);
        
        res = norm(b - A*x);
        resvec = [resvec,res];
        
    end
    t2 = toc();
    if verbose
        fprintf("system size %d \n", size(A,1));
        fprintf("mesh subdivisions %d \n", n + 1);
        fprintf("solving exactly at level %d with mesh subd %d \n", max_lvl, (n+1)/(2^(max_lvl - 1)));
        fprintf("nu1=%d nu2=%d \n", nu1,nu2);
        fprintf("time to prepare operators and grids %e \n", t1)
        fprintf("time to solution %e \n",t2);
        
        
    end
    
end

