function  P = get_PR_opt(n_2h)
%coarse to finer grid interpolation procedure
%uses linear interpolation

   n_h = 2*n_2h + 1;
   %P = spalloc(n_h, n_2h, n_2h * 3);
   row = zeros(3*n_2h,1);
   col = zeros(3*n_2h,1);
   val = zeros(3*n_2h,1);
   i = 1;
   for i_2h = 1:(n_2h)
        i_h = 2*i_2h;
        
        row(i:i+2) = [i_h - 1, i_h, i_h + 1];
        col(i:i+2) = [i_2h, i_2h, i_2h];
        val(i:i+2) = [1,2,1];
        
        
        %P(i_h - 1, i_2h) = 1;
        %P(i_h + 1, i_2h) = 1;
        %P(i_h, i_2h) = 2;
        
        i = i + 3;
               
   end     
   
   P = sparse(row,col, val, n_h, n_2h);
   
end

