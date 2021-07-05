function  [P,R] = get_PR_fw(nx_2h)
%interpolation function, works with a 2d grid, uses linear interpolation
%to make the prolongation into a finer grid with 2h as grid spacing

   nx_h = nx_2h*2 + 1;
   
   P = spalloc(nx_h^2, nx_2h^2, nx_2h^2 * 9) ;
   
   for i = 1:nx_2h
       ii = 2*i;
        for j = 1:nx_2h
            jj = 2*j;
            idx_2h = (i - 1)*nx_2h + j; 
            idx_h = (ii - 1)*nx_h + jj;
            
            P(idx_h,idx_2h) = 4;
            
            P(idx_h - nx_h - 1, idx_2h) = 1;
            P(idx_h + nx_h - 1, idx_2h) = 1;
            P(idx_h - nx_h + 1, idx_2h) = 1;
            P(idx_h + nx_h + 1, idx_2h) = 1;
            
            P(idx_h - nx_h , idx_2h) = 2;
            P(idx_h + nx_h , idx_2h) = 2;
            P(idx_h - 1, idx_2h) = 2;
            P(idx_h  + 1, idx_2h) = 2;
            
            
        end
   end
   
   R = P.';
   
   R = R/16;
   
   P = P/4;
   
   
   

end

