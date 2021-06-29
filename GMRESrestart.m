function [x_k,iter,resvec,flag] = GMRESrestart(A,b,outer,inner,x0,tol)
flag=0;
resvec=[];
b2=tol*b;


for k=1:outer
    [x_k,iter_k,resvec_k,flag] = mygmres(A,b,tol,inner,x0);
    rk=resvec_k(end);
    resvec=[resvec;rk];
    if ( rk < b2 )
        flag=0; %convergenza
        return
    else
        x0=x_k;
        flag = 1; %nonconvergenza
    end
end
iter=(k - 1)*inner + iter_k;


end







% function [x,iter,resvec,flag] = restartgmresVERO(A,b,tol,p,maxit,x0)
% 
% flag=0;
% nb=norm(b);
% if nb>0 
%     b2=tol*nb;
% else
%     b2=tol;
% end
% n = size(A,1);
% resvec=[];
% V = zeros(n,p+1); 
% H = zeros(p+1,n);
% iter=0;
% 
% 
% for iter=1:maxit
%     
%     r=b-A*x0;
%     rho=norm(r); 
%     beta=rho;
%     V(:,1) = r/beta;
%     iter=iter+1;
%     
%     for i=1:p
%         V(:,i+1)=A*V(:,i);
%         iter=iter+1;
%         for j=1:i
%             H(j,i)=V(:,i+1)'*V(:,j);
%             V(:,i+1)=V(:,i+1)-H(j,i)*V(:,j);
%         end
%         H(i+1,iter)=norm(V(:,i+1));
%         V(:,i+1)=V(:,i+1)/H(i+1,i);
% 
%         [Q,R] = qr(H(1:i+1,1:i));
%         if H(i+1,i) == 0, break, end
%         rho=abs(beta*Q(1,i+1)); 
%         resvec=[resvec;rho];
% 
%     end
%     y = zeros(i+1,1);
%     y(1:i,1) = R(1:i,:)\(beta*Q(1,1:i)');
%     x = x0+V(:,1:i+1)*y;
%     %error = norm(x-x_exact);
%     r=norm(b-A*x);
%     if ( r <= b2 )
%         return
%     else
%         x0=x;
%         flag = 0; 
%     end
% end
% end
%         


