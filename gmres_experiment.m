function gmres_experiment()
    matrices = {"olm100", "oscil_dcop_02", "cavity05", "coater2", "Dubcova1"};

    ks_per_mat = {
        [5, 10, 20, 30],
        [10, 25, 50, 100],
        [100, 200, 350, 500],
        [250, 500, 650, 800],
        [250, 500, 650, 800]
    };

    tols = [10e-6, 10e-7, 10e-8, 10e-9, 10e-10, 10e-11];

    maxit_per_mat = {
        [100, 200, 300, 400, 500, 600, 700],
        [100, 200, 300, 400, 500, 600, 700],
        [100, 200, 300, 400, 500, 600, 700, 800],
        [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000],
        [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    };

    for i = 1 : numel(matrices)
        matrices{i}
        load(sprintf("in/%s.mat", matrices{i}));
        A = Problem.A;
        [x, flag, relres, iter, resvec, tol, maxit] = test_gmres_params(A, ones(rows(A), 1), ks_per_mat{i}, tols, maxit_per_mat{i});

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
