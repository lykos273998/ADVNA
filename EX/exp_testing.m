%testing multigrid solver against cg 
%rhs s.t. x_exc is ones
%x0 
addpath("1D","2D");
lvl_1d = [2 3 4];

test_1d = 1;
test_2d = 1;

meshes_1d = [32 64 128 256 512 1024];

lvl_2d = [2 3];

meshes_2d = [32 64 128 256, 512];

outfile = fopen('results.txt','w');



maxit = 2000;
tol = 1e-10;
if test_1d 
    
    fprintf(outfile,"*******************************\n");
    fprintf(outfile,"*            1D mesh          *\n");
    fprintf(outfile,"*******************************\n\n");
    
    
    for max_lvl = lvl_1d
        fprintf(outfile,"\n++++++++++++++++++++++\n");
        fprintf(outfile,"V-Cycle levels %d \n\n", max_lvl);
        for N = meshes_1d
            fprintf(outfile,"\t Mesh subdivisions (1D) %d \n",N);
            fprintf(outfile,"\t -------------------------- \n");
            A = L1D(N);
            x_exc = ones(size(A,1),1);
            x0 = zeros(size(A,1),1);
            b = A*x_exc;
            tic();
            [xmg, img, resmg ] = MG_1D_gen(A,b,x0,maxit,tol,max_lvl);
            tmg = toc();

            tic();
            [xcg, icg, rescg] = mypcg(A,b,x0,maxit,tol,speye(size(A,1)));
            tcg = toc();

            rate_mg = -log10(norm(x_exc - xmg)/norm(x_exc - x0))/img;
            rate_cg = -log10(norm(x_exc - xcg)/norm(x_exc - x0))/icg;

            fprintf(outfile, "\t\t MultiGrid >>> Iter: %d  Conv Rate: %.4f  err norm: %e CPU: %e\n", img, rate_mg, norm(x_exc - xmg), tmg);
            fprintf(outfile, "\t\t CG        >>> Iter: %d  Conv Rate: %.4f  err norm: %e CPU: %e\n\n", icg, rate_cg, norm(x_exc - xcg), tcg);





        end

    end
end


if test_2d
    fprintf(outfile,"\n\n\n");
    fprintf(outfile,"*******************************\n");
    fprintf(outfile,"*            2D mesh          *\n");
    fprintf(outfile,"*******************************\n\n");

    maxit = 2000;
    tol = 1e-10;

    for max_lvl = lvl_2d
        fprintf(outfile,"\n++++++++++++++++++++++\n");
        fprintf(outfile,"V-Cycle levels %d \n\n", max_lvl);
        for N = meshes_2d
            fprintf(outfile,"\t Mesh subdivisions (2D) %d \n",N);
            fprintf(outfile,"\t -------------------------- \n");
            
            A = delsq(numgrid('S',N + 1));
            x_exc = ones(size(A,1),1);
            x0 = zeros(size(A,1),1);
            b = A*x_exc;
            tic();
            [xmg, img, resmg ] = MG_2D_gen(A,b,x0,maxit,tol,N-1,max_lvl);
            tmg = toc();

            tic();
            [xcg, icg, rescg] = mypcg(A,b,x0,maxit,tol,speye(size(A,1)));
            tcg = toc();

            rate_mg = -log10(norm(x_exc - xmg)/norm(x_exc - x0))/img;
            rate_cg = -log10(norm(x_exc - xcg)/norm(x_exc - x0))/icg;

            fprintf(outfile, "\t\t MultiGrid >>> Iter: %d  Conv Rate: %.4f  err norm: %e CPU: %e\n", img, rate_mg, norm(x_exc - xmg), tmg);
            fprintf(outfile, "\t\t CG        >>> Iter: %d  Conv Rate: %.4f  err norm: %e CPU: %e\n\n", icg, rate_cg, norm(x_exc - xcg), tcg);





        end

    end
end


