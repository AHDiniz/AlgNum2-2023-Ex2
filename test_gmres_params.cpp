#include "octave/oct.h"
#include "octave/parse.h"
#include <vector>
#include <limits>

void gmres(Matrix &A, ColumnVector &b, int k, double tol, int max_iter, octave_value_list &result)
{
    octave_value_list gmresin;
    gmresin(0) = A;
    gmresin(1) = b;
    gmresin(2) = k;
    gmresin(3) = tol;
    gmresin(4) = max_iter;

    result = octave::feval("gmres", gmresin);
}

struct GMRESInput
{
    int k;
    double tol;
    int max_iter;
};

DEFUN_DLD(test_gmres_params, args, nargout, "Tests the best parameters for the GMRES method.")
{
    Matrix A = args(0).matrix_value();
    ColumnVector b = args(1).column_vector_value();
    ColumnVector ks = args(2).column_vector_value();
    ColumnVector tols = args(3).column_vector_value();
    ColumnVector maxits = args(4).column_vector_value();

    std::vector<GMRESInput> inputs;

    for (int i = 0; i < ks.numel(); ++i)
    {
        for (int j = 0; j < tols.numel(); ++j)
        {
            GMRESInput in;
            in.k = ks(i);
            in.tol = tols(j);
            in.max_iter = tols(l);
        }
    }

    std::vector<octave_value_list> results;
    for (GMRESInput in : inputs)
    {
        octave_value_list result;
        gmres(A, b, in.k, in.tol, in.max_iter, result);
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
            retval(7) = inputs[i].tol;
            retval(8) = inputs[i].max_iter;
        }
    }

    return retval;
}
