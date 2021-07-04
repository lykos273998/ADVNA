function  r_h = c_to_f(r_2h)
%coarse to finer grid interpolation procedure
%uses linear interpolation

   n = size(r_2h,1);
   r_h = zeros(n*2+1,1);
   j = 2*1;
    r_h(j- 1) = 0.5 *( r_2h(1));
    r_h(j) = r_2h(1);
    r_h(j+1) = 0.5 *( r_2h(1) + r_2h(2) );
    
   for i = 2:(n-1)
        j = 2*i;
        r_h(j- 1) = 0.5 *( r_2h(i) + r_2h(i -1));
        r_h(j) = r_2h(i);
        r_h(j+1) = 0.5 *( r_2h(i) + r_2h(i + 1));
        
   end     
   j = 2*n;
    r_h(j - 1) = 0.5 * ( r_2h(n) + r_2h(n-1) );
    r_h(j) = r_2h(n);
    r_h(j+1) = 0.5 *( r_2h(n) );
end

