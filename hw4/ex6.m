%bratu problem setup
h = 5e-3;
N = 1/h;
B = (-1)*delsq(numgrid('S',N+1));

n = size(B,1);
lambda = 6.5;

F = @(u) B*u + h^2 * lambda * exp(u);
Jac = @(u) B + spdiags(h^2 * lambda * exp(u),0,n,n);

tol = 1e-13;
x0 = 0.1*ones(n,1);

outfile = fopen("ex6_results.txt",'w');
itmax = 100;
%newton w direct solution
fprintf(outfile,"/*************************************\n");
fprintf(outfile," * Classical Newton                  *\n");
fprintf(outfile," *************************************/\n");

tic();
[x_newt, iter_newt, resvec_newt] = newton(x0,F,Jac, tol, itmax,1);
t = toc();

fprintf(outfile, "iteration\t||F(x_k)||\n");
for i=1:iter_newt
    fprintf(outfile,"   %d    \t %e \n", i, resvec_newt(i));
end

fprintf(outfile,"CPU TIME: %.2f\n", t);
figure();
semilogy(resvec_newt,'.-');
hold on
fprintf(outfile,"--------------------------------------\n\n");

fprintf(outfile,"/*************************************\n");
fprintf(outfile," * Quasi Newton                      *\n");
fprintf(outfile," *************************************/\n");

J0 = Jac(x0);
B0 = @(x) J0;
tic();

[x_quasi, iter_quasi, resvec_quasi] = quasi_newton(x0,F,B0, tol, itmax);
t = toc();
fprintf(outfile, "iteration\t||F(x_k)||\n");
for i=1:iter_quasi
    fprintf(outfile,"   %d    \t %e \n", i, resvec_quasi(i));
end
fprintf(outfile,"CPU TIME: %.2f\n", t);
semilogy(resvec_quasi,'.-')
legend("classical newton","quasi newton");
xlabel("iteration");
ylabel("||F(x_k)||");
hold off
fprintf(outfile,"--------------------------------------\n\n");


fprintf(outfile,"/*************************************\n");
fprintf(outfile," * Inexact Newton: case i            *\n");
fprintf(outfile," *************************************/\n");

eta_max = 0.1;

exit_tol = tol*norm(F(x0));
res = exit_tol + 1;
resvec = [];
iter_inx = 0;
x_inx = x0;
f = -F(x0);

gmres_restart = 50;
gmres_maxit = 10;
setup.type = 'ilutp';
setup.droptol = 1e-2;
fprintf(outfile, "iteration	||F(x_k)||	     eta 	       Lin IT\n");
tic();
while res > exit_tol & iter_inx < itmax
        iter_inx = iter_inx + 1;
        J = Jac(x_inx);
        [L,U] = ilu(J,setup);
        eta = eta_max;
        [s,a,b,it_gmres] = gmres(J, f, gmres_restart, eta, gmres_maxit, L, U, x_inx);
        x_inx = x_inx + s;

        f = -F(x_inx);                
        res = norm(f);

        resvec = [resvec,res];
        fprintf(outfile,"   %d    \t %e \t %e \t %d \n", iter_inx, res, eta, it_gmres(1)*gmres_restart + it_gmres(2));
end
t = toc();
fprintf(outfile,"CPU TIME: %.2f\n", t);
figure();
semilogy(resvec,'.-');
hold on
fprintf(outfile,"--------------------------------------\n\n");


fprintf(outfile,"/*************************************\n");
fprintf(outfile," * Inexact Newton: case ii           *\n");
fprintf(outfile," *************************************/\n");

eta_max = 0.1;

exit_tol = tol*norm(F(x0));
res = exit_tol + 1;
resvec = [];
iter_inx = 0;
x_inx = x0;
f = -F(x0);

fprintf(outfile, "iteration	||F(x_k)||	     eta 	       Lin IT\n");
eta = eta_max;

tic();
while res > exit_tol & iter_inx < itmax
        iter_inx = iter_inx + 1;
        J = Jac(x_inx);
        [L,U] = ilu(J,setup);
        eta = eta/3;
        [s, a,b,it_gmres] = gmres(J, f, gmres_restart, eta, gmres_maxit, L, U, x_inx);
        x_inx = x_inx + s;

        f = -F(x_inx);                
        res = norm(f);

        resvec = [resvec,res];
        fprintf(outfile,"   %d    \t %e \t %e \t %d \n", iter_inx, res, eta, it_gmres(1)*gmres_restart + it_gmres(2));
end
t = toc();

fprintf(outfile,"CPU TIME: %.2f\n", t);
semilogy(resvec,'.-');
hold on
fprintf(outfile,"--------------------------------------\n\n");



fprintf(outfile,"/*************************************\n");
fprintf(outfile," * Inexact Newton: case iii          *\n");
fprintf(outfile," *************************************/\n");

l = 0.95;
eta_max = 0.1;

exit_tol = tol*norm(F(x0));
res = exit_tol + 1;
resvec = [];
iter_inx = 0;
x_inx = x0;
f = -F(x0);

fprintf(outfile, "iteration	||F(x_k)||	     eta 	       Lin IT\n");

tic();
eta = eta_max;

%first iteration
iter_inx = iter_inx + 1;
J = Jac(x_inx);
[L,U] = ilu(J,setup);
[s, a,b,it_gmres] = gmres(J, f, gmres_restart, eta, gmres_maxit, L, U, x_inx);
x_inx = x_inx + s;

f = -F(x_inx);                
res = norm(f);

resvec = [resvec,res];
fprintf(outfile,"   %d    \t %e \t %e \t %d \n", iter_inx, res, eta, it_gmres(1)*gmres_restart + it_gmres(2));


while res > exit_tol & iter_inx < itmax
        iter_inx = iter_inx + 1;
        J = Jac(x_inx);
        [L,U] = ilu(J,setup);
        eta = min(eta_max, l*res);
        [s, a,b,it_gmres] = gmres(J, f, gmres_restart, eta, gmres_maxit, L, U, x_inx);
        x_inx = x_inx + s;

        f = -F(x_inx);                
        res = norm(f);

        resvec = [resvec,res];
        fprintf(outfile,"   %d    \t %e \t %e \t %d \n", iter_inx, res, eta, it_gmres(1)*gmres_restart + it_gmres(2));
end
t = toc();


fprintf(outfile,"CPU TIME: %.2f\n", t);
semilogy(resvec,'.-');
hold on
fprintf(outfile,"--------------------------------------\n\n");



fprintf(outfile,"/*************************************\n");
fprintf(outfile," * Inexact Newton: case  iv          *\n");
fprintf(outfile," *************************************/\n");

eta_max = 0.1;

exit_tol = tol*norm(F(x0));
res = exit_tol + 1;
resvec = [];
iter_inx = 0;
x_inx = x0;
f = -F(x0);

fprintf(outfile, "iteration	||F(x_k)||	     eta 	       Lin IT\n");
tic();
eta = eta_max;


while res > exit_tol & iter_inx < itmax
        iter_inx = iter_inx + 1;
        J = Jac(x_inx);
        [L,U] = ilu(J,setup);
        
        if (iter_inx > 2)
            eta = min(eta_max, l*(resvec(end)/resvec(end-1))^2);
        else
            eta = eta_max;
        end
        
        [s, a,b,it_gmres] = gmres(J, f, gmres_restart, eta, gmres_maxit, L, U, x_inx);
        x_inx = x_inx + s;

        f = -F(x_inx);                
        res = norm(f);

        resvec = [resvec,res];
        fprintf(outfile,"   %d    \t %e \t %e \t %d \n", iter_inx, res, eta, it_gmres(1)*gmres_restart + it_gmres(2));
end

t = toc();

fprintf(outfile,"CPU TIME: %.2f\n", t);
semilogy(resvec,'.-');
hold on
semilogy(resvec_newt,'.-');
legend("inx newton: case i","inx newton: case ii","inx newton: case iii","inx newton: case iv","classical newton");
xlabel("iteration");
ylabel("||F(x_k)||");
hold off

fprintf(outfile,"--------------------------------------\n\n");
