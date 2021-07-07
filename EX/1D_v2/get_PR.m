function  P = get_PR(n_2h)
%coarse to finer grid interpolation procedure
%uses linear interpolation

   n_h = 2*n_2h + 1;
   P = spalloc(n_h, n_2h, n_2h * 3);
    
   for i_2h = 1:(n_2h)
        i_h = 2*i_2h;
        P(i_h - 1, i_2h) = 1;
        P(i_h + 1, i_2h) = 1;
        P(i_h, i_2h) = 2;
               
   end     
   
end

