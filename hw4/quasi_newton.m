function [x,iter,resvec] = quasi_newton(x0,F,B0,tol, itmax)
%QUASI_NEWTON Summary of this function goes here
%   Detailed explanation goes here
    
    f = -F(x0);
    [L,U] = lu(B0(x0));
    s = U\(L\f);
    %s = B0(x0)\f;
    k = 0;
    exit_tol = tol*norm(F(x0));
    err = norm(f);
    resvec = [];
    s_mat = [];
    x = x0; 
    
    
    while k < itmax && err > exit_tol
        %s_mat
        x = x + s;
        s_mat = [s_mat, s];
        f = -F(x);
        [L,U] = lu(B0(x));
        z = U\(L\f);
        
        %z = B0(x)\f;
                    
        for j = 1:size(s_mat,2)-1
            %disp(k)
            %disp(size(s_mat,2))
            z = z + s_mat(:,j+1) * (s_mat(:,j).'*z)/(norm(s_mat(:,j))^2);
        end
        
        s = z/(1 - s.'*z/(norm(s)^2));
        
            
        k = k + 1;
        
        err = norm(f);
        resvec = [resvec, err];
    end
    
    iter = k;
    
end

