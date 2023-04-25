function pcg_experiment()
    # matrices = {"mesh3em5", "plat362", "662_bus", "s1rmq4m1", "bcsstk36", "pdb1HYS", "Dubcova3"};
    # matrices = {"mesh3em5", "plat362", "662_bus"}
    matrices = {"mesh3em5", "plat362", "662_bus", "s1rmq4m1", "bcsstk36"};

    tols = [10e-6, 10e-7, 10e-8, 10e-9, 10e-10, 10e-11];

    maxits = [1000, 5000, 7500, 10000];

    for i = 1 : numel(matrices)
        matrices{i}
        load(sprintf("in/%s.mat", matrices{i}));
        A = Problem.A;
        [x, flag, relres, iter, resvec, tol, maxit] = test_pcg_params(A, A * ones(rows(A), 1), tols, maxits);

        f = fopen(sprintf("out/%s.txt", matrices{i}), "w");
        fprintf(f, "n = %d\n", rows(A));
        fprintf(f, "nnz = %d\n", nnz(A));
        fprintf(f, "flag = %d\n", flag);
        fprintf(f, "iterations = %d\n", iter);
        fprintf(f, "sol. norm = %f\n", norm(x, inf));
        fprintf(f, "tol = %f\n", tol);
        fprintf(f, "maxit = %d\n", maxit);
        fclose(f);

        hf = figure();
        plot(1:numel(resvec), resvec);
        print(hf, sprintf("out/%s.png", matrices{i}), "-dpng");
    endfor
end
