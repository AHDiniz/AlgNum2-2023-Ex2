#include "octave/oct.h"
#include "octave/parse.h"
#include <vector>
#include <limits>

void pcg(Matrix &A, ColumnVector &b, double tol, int max_iter, octave_value_list &result)
{
    octave_value_list cgin;
    cgin(0) = A;
    cgin(1) = b;
    cgin(2) = tol;
    cgin(3) = max_iter;

    result = octave::feval("pcg", cgin);
}

struct CGInput
{
    double tol;
    int max_iter;
};

DEFUN_DLD(test_pcg_params, args, nargout, "Tests the best parameters for the Conjugate Gradient method.")
{
    Matrix A = args(0).matrix_value();
    ColumnVector b = args(1).column_vector_value();
    ColumnVector tols = args(2).column_vector_value();
    ColumnVector maxits = args(3).column_vector_value();

    std::vector<CGInput> inputs;

    for (int i = 0; i < tols.numel(); ++i)
    {
        for (int j = 0; j < maxits.numel(); ++j)
        {
            CGInput cgin;
            cgin.tol = tols(i);
            cgin.max_iter = maxits(j);
            inputs.push_back(cgin);
        }
    }

    /*unsigned max_threads = std::thread::hardware_concurrency();
    max_threads = max_threads / 2 + 1;
    int executions = 0;
    std::vector<std::thread> threads;
    std::vector<octave_value_list> results;
    while (executions < inputs.size())
    {
        octave_value_list result;
        std::thread t(&pcg, A, b, inputs[executions].tol, inputs[executions].max_iter, std::ref(result));
        results.push_back(result);
        threads.push_back(t);
        if (threads.size() == max_threads)
        {
            for (int i = 0; i < threads.size(); ++i)
            {
                threads[i].join();
                ++executions;
            }
            threads.clear();
        }
    }*/
    std::vector<octave_value_list> results;
    for (CGInput cgin : inputs)
    {
        octave_value_list result;
        pcg(A, b, cgin.tol, cgin.max_iter, result);
        results.push_back(result);
    }

    double best_relres = std::numeric_limits<double>::max();
    int best_flag = 4;
    octave_value_list retval;
    for (int i = 0; i < results.size(); ++i)
    {
        octave_value_list result = results[i];
        int flag = result(1).int_value();
        int relres = result(2).double_value();

        if (flag < best_flag && relres < best_relres)
        {
            retval = result;
            retval(6) = inputs[i].tol;
            retval(7) = inputs[i].max_iter;
        }
    }

    return retval;
}
