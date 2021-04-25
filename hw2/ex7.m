NS = [201, 401];

tol = 1e-8;
maxiter = 5000;
p = 10;
for N = NS
    A = delsq(numgrid('S', N + 1));
    n = size(A,1);
    x_exact = ones(n,1);
    b = A*x_exact;
    
    L = ichol(A);
    
    [W, Lambda] = eigs(A, L*L', p, 'sm');
    H = sparse(W' * A* W);
    H = inv(H);
    
    P_handle = @(x) L'\(L\x) + W*(H*(W'*x));
    
    [x1, flag1, relres1, iter1, resvec1] = pcg(A,b,tol,maxiter, L, L');
    [x2, flag2, relres2, iter2, resvec2] = pcg(A,b,tol,maxiter, P_handle);
    figure()
    semilogy(1:(iter1+1), resvec1,'*-', 1:(iter2+1), resvec2, '*-')
    xlabel('Iterations');
    ylabel("||r_k||");
    legend("IC(0)", "Spectral preconditioner");
    title(sprintf("N = %d, System size = %d",N,n));
    
end


