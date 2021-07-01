function e_h = V_cycle(r_h, nx, nu1, nu2, A_vec, tol, lvl, max_lvl)
    if lvl == max_lvl
        e_h = A_vec{lvl}\r_h;
    else
        e_h_0 = zeros(size(r_h));
        e_h = SOR(A_vec{lvl},r_h,e_h_0,nu1,tol,1);
        
        nx_2 = floor(nx/2);
        r_2h = f_to_c_2d(r_h - A_vec{lvl}*e_h,nx_2);
        e_2h = V_cycle(r_2h,nx_2,nu1,nu2,A_vec,tol,lvl+1,max_lvl);
        
        e_h = e_h + c_to_f_2d(e_2h,nx_2);
        e_h = SOR(A_vec{lvl},r_h,e_h,nu2,tol,1);
                
    end
    
end