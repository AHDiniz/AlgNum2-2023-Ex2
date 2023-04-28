function gmres_experiment()
    % matrices = {"olm100", "oscil_dcop_02", "cavity05", "coater2", "Dubcova1"};
    matrices = {"coater2", "Dubcova1"};

    ks = [5, 25, 50, 75];

    % tols = [10e-6, 10e-8, 10e-9, 10e-11];
    tols = [10e-8];

    for i = 1 : numel(matrices)
        matrices{i}
        load(sprintf("in/%s.mat", matrices{i}));
        A = Problem.A;
        [x, cflag, relres, iter, resvec, k, tol, maxit] = test_gmres_params(A, ones(rows(A), 1), ks, tols, 1e5);

        n_iter = iter(1,1) * k + iter(1,2);

        f = fopen(sprintf("out/%s.txt", matrices{i}), "w");
        fprintf(f, "n = %d\n", rows(A));
        fprintf(f, "nnz = %d\n", nnz(A));
        fprintf(f, "flag = %d\n", cflag);
        fprintf(f, "iterations = %d\n", n_iter);
        fprintf(f, "sol. norm = %f\n", norm(x, inf));
        fprintf(f, "k = %d\n", k);
        fprintf(f, "tol = %f\n", tol);
        fprintf(f, "maxit = %d\n", maxit);
        fclose(f);

        m = min(n_iter, numel(resvec));
        aux_resvec = zeros(m, 1);
        for j = 1 : m
            aux_resvec(j) = resvec(j);
        end

        hf = figure();
        plot(1:m, log(aux_resvec));
        print(hf, sprintf("out/%s.png", matrices{i}), "-dpng");
    endfor
end
