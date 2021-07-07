function  [P] = get_PR_fw(nx_2h)
%interpolation function, works with a 2d grid, uses linear interpolation
%to make the prolongation into a finer grid with 2h as grid spacing

   nx_h = nx_2h*2 + 1;
   
   %P = spalloc(nx_h^2, nx_2h^2, nx_2h^2 * 9) ;
   
   row = zeros(nx_2h^2 * 9,1);
   col = zeros(nx_2h^2 * 9,1);
   val = zeros(nx_2h^2 * 9,1);
   mat_idx = 1;
   for i = 1:nx_2h
       ii = 2*i;
        for j = 1:nx_2h
            jj = 2*j;
            idx_2h = (i - 1)*nx_2h + j; 
            idx_h = (ii - 1)*nx_h + jj;
            row(mat_idx:mat_idx+8) = [idx_h, idx_h - nx_h - 1, idx_h + nx_h - 1, idx_h - nx_h + 1, idx_h + nx_h + 1, idx_h - nx_h , idx_h + nx_h, idx_h - 1, idx_h + 1];
            col(mat_idx:mat_idx+8) = [idx_2h, idx_2h,idx_2h,idx_2h,idx_2h,idx_2h,idx_2h,idx_2h,idx_2h,];
            val(mat_idx:mat_idx+8) = [4 1 1 1 1 2 2 2 2];
            
            mat_idx = mat_idx+9;
            %P(idx_h,idx_2h) = 4;
            
            %P(idx_h - nx_h - 1, idx_2h) = 1;
            %P(idx_h + nx_h - 1, idx_2h) = 1;
            %P(idx_h - nx_h + 1, idx_2h) = 1;
            %P(idx_h + nx_h + 1, idx_2h) = 1;
            
            %P(idx_h - nx_h , idx_2h) = 2;
            %P(idx_h + nx_h , idx_2h) = 2;
            %P(idx_h - 1, idx_2h) = 2;
            %P(idx_h  + 1, idx_2h) = 2;
            
            
        end
   end
   
   
   
   P = sparse(row,col,val,nx_h^2,nx_2h^2);
   
   
   

end

