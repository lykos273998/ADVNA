function [x, iter,resvec] = MG_2D(A,b,x0,maxit,tol,nx)
    %2 level multigrid test
    iter = 0;
    exit_tol = tol*norm(b);
    x = x0;
    res = norm(b - A*x);
    resvec = [ res ];
    
    nx2 = floor(nx/2);
    nx4 = floor(nx2/2);
    
    n = size(A,1);
    A2 = delsq(numgrid('S',nx2+2));
    A4 = delsq(numgrid('S',nx4+2));
    
    n2 = size(A2,1);
    n4 = size(A4,1);
    
    
   
    
    nu1 = 5;
    nu2 = 3;
    
    r_2h = zeros(n2,1);
    e_2h_0 = zeros(n2,1);
    e_4h_0 = zeros(n4,1);
    nu_tot = 0;
    
    while iter < maxit && res > exit_tol
        iter = iter + 1;
        
        %first level
        [x, i, vdiff] = SOR(A,b,x,nu1,tol,1);
        
        r_h = b - A*x;
      
        r_2h = f_to_c_2d(r_h,nx2);
        
        %second level        
        [e_2h, i, vdiff] = SOR(A2,r_2h,e_2h_0,nu1,tol,1);
        
        
        r_4h = f_to_c_2d(r_2h - A2*e_2h,nx4);
        
        %third level
        %[e_4h, i, vdiff] = SOR(A4,r_4h,e_4h_0,nu2,tol,1);
        e_4h = A4\r_4h;
        
        e_2h = e_2h + c_to_f_2d(e_4h,nx4);
        
        [e_2h, i, vdiff] = SOR(A2,r_2h,e_2h,nu2,tol,1);
      
        e_h = c_to_f_2d(e_2h,nx2);
        
        x = x + e_h;
        
      
        [x, i, vdiff] = SOR(A,b,x,nu2,tol,1);
        
        res = norm(b - A*x);
        resvec = [resvec,res];
        
        
    end
    
end

