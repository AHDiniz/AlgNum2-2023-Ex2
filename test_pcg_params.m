function [best_x, best_flag, best_relres, best_iter, best_resvec, best_tol, best_maxit] = test_pcg_params(A, b, tols, maxits)

    inputs = cell(numel(tols) * numel(maxits), 2);
    for i = 1 : numel(tols)
        for j = 1 : numel(maxits)
            index = (i - 1) * numel(maxits) + j;
            inputs{index,1} = tols(i);
            inputs{index,2} = maxits(j);
        endfor
    endfor

    best_sol_diff = inf;
    for i = 1 : numel(tols) * numel(maxits)
        tol = inputs{i,1};
        maxit = inputs{i,2};
        [x, cflag, relres, iter, resvec] = conjgrad(A, b, tol, maxit);
        sol_diff = abs(abs(norm(x, inf)) - 1);
        if sol_diff < best_sol_diff
            best_x = x;
            best_flag = cflag;
            best_relres = relres;
            best_iter = iter;
            best_resvec = resvec;
            best_tol = tol;
            best_maxit = maxit;
            best_sol_diff = sol_diff;
        endif
    endfor
end
