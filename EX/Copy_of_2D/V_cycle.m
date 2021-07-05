function e_h = V_cycle(r_h, nu1, nu2, A_vec,P_vec,R_vec, tol, lvl, max_lvl)
%V-Cycle function,
%calls recursively itself up to the
%maximium depth specified,once reached the coarser grid it solves exactly
%the equation Ae = r
    if lvl == max_lvl
        e_h = A_vec{lvl}\r_h;
    else
        e_h_0 = zeros(size(r_h));
        e_h = SOR(A_vec{lvl},r_h,e_h_0,nu1,tol,1);
        
        r_2h = R_vec{lvl+1} * (r_h - A_vec{lvl}*e_h);
        e_2h = V_cycle(r_2h, nu1,nu2,A_vec,P_vec,R_vec,tol,lvl+1,max_lvl);
        
        e_h = e_h + P_vec{lvl + 1} * e_2h;
        e_h = SOR(A_vec{lvl},r_h,e_h,nu2,tol,1);
                
    end
    
end