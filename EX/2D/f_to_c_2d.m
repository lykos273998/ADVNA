function r_2h = f_to_c_2d(r_h,nx)
%restriction function
%injects values of a finer grid into a coarser one,
%uses simple injection of pouints (2i,2j) into (i,j) of the finer grid
   
   r_2h = zeros(nx^2,1);
  
   nx_h = 2 * nx + 1;
      
   for i = 1:nx
       ii = 2*i;
        for j = 1:nx
                        
            r_2h((i - 1)*nx + j) = r_h( (ii - 1)*nx_h + jj);
            
        end
   end
   
   r_2h = r_2h;
end

