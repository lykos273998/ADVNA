function  r_h = c_to_f_2d(r_2h,nx)
   n = size(r_2h,1);
   r_h = zeros((nx*2+1)^2,1);
    
   nx_h = 2*nx + 1;
   
   for i = 1:nx
       ii = 2*i;
        for j = 1:nx
            jj = 2*j;
            val = r_2h((i - 1)*nx + j);             
            r_h( (ii - 1)*nx_h + jj) = r_h( (ii - 1) *nx_h + jj) + 4*val;
            
            r_h( (ii - 1 - 1)*nx_h + (jj - 1)) = r_h( (ii - 1 - 1)*nx_h + (jj - 1)) + val;
            r_h( (ii + 1 - 1)*nx_h + (jj - 1)) = r_h( (ii + 1 - 1)*nx_h + (jj - 1)) + val;
            r_h( (ii - 1 - 1)*nx_h + (jj + 1)) = r_h( (ii - 1 - 1)*nx_h + (jj + 1)) + val;
            r_h( (ii + 1 - 1)*nx_h + (jj + 1)) = r_h( (ii + 1 - 1)*nx_h + (jj + 1)) + val;
            
            r_h( (ii - 1 )*nx_h + (jj - 1)) = r_h( (ii - 1)*nx_h + (jj - 1)) + 2*val;
            r_h( (ii - 1)*nx_h + (jj + 1)) = r_h( (ii - 1)*nx_h + (jj + 1)) + 2*val;
            r_h( (ii - 1 - 1)*nx_h + (jj)) = r_h( (ii - 1 - 1)*nx_h + (jj)) + 2*val;
            r_h( (ii + 1 - 1)*nx_h + (jj)) = r_h( (ii + 1 - 1)*nx_h + (jj)) + 2*val;
            
        end
   end
   
   r_h = r_h/4;
   
end

