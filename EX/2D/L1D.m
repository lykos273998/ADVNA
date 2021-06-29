function L = L1D(N)
    u= -2*ones(N-1,1);
    w=ones(N-2,1);
    L= sparse(diag(u)+diag(w,-1)+diag(w,1));
end

