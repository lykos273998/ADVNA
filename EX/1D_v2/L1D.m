function L = L1D(N)
%returns 1d discrete laplacian
    u= -2*ones(N-1,1);
    w=ones(N-2,1);
    L= sparse(diag(u)+diag(w,-1)+diag(w,1));
end

