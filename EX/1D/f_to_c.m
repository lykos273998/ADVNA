function r_2h = f_to_c(r_h)
    n = floor(size(r_h,1)/2);
    r_2h = zeros(n,1);
    for i = 1:n
        j = 2*i;
        r_2h(i) = 0.25*(r_h(j-1) + 2*r_h(j) + r_h(j+1));
    end  
    
end

