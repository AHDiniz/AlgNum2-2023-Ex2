#include "octave/oct.h"
#include "octave/parse.h"
#include <vector>
#include <cmath>

inline double dot_product(const ColumnVector &a, const ColumnVector &b)
{
    double d = 0.0;
    for (int i = 0; i < a.numel(); ++i)
        d += a(i) * b(i);
    return d;
}

inline double euc_norm(const ColumnVector &a)
{
    return sqrt(dot_product(a, a));
}

DEFUN_DLD(conjgrad, args, nargout, "C++ implementation of the Conjugate Gradient method.")
{
    Matrix A = args(0).matrix_value();
    ColumnVector b = args(1).column_vector_value();
    double tol = args(2).double_value();
    int max_iter = args(3).int_value();

    const double EPS = octave::feval("eps")(0).double_value();

    int n = A.rows();

    std::vector<double> residuals;

    ColumnVector x0(n);
    for (int i = 0; i < n; ++i)
        x0(i) = 0.0;
    ColumnVector d0 = b;
    ColumnVector r0 = b;
    ColumnVector d(n);
    ColumnVector x(n);
    ColumnVector r(n);

    double delta = dot_product(r0, r0);
    double delta0 = delta;

    int iter = 0;
    int flag = 0;

    while (delta > tol * tol * delta0 && iter < max_iter)
    {
        ColumnVector v = A * d0;
        double lambda = delta / dot_product(d0, v);
        x = x0 + lambda * d0;
        r = r0 - lambda * v;
        residuals.push_back(euc_norm(r));
        delta0 = delta;
        delta = dot_product(r, r);
        double beta = delta0 / delta;
        d = r + beta * d0;
        ++iter;

        if (euc_norm(x - x0) <= EPS * euc_norm(x))
        {
            flag = 3;
            break;
        }

        x0 = x;
        r0 = r;
        d0 = d;
    }

    flag = (flag != 3 && iter >= max_iter) ? 1 : flag;
    double relres = residuals[residuals.size() - 1] / euc_norm(b);
    ColumnVector resvec(residuals.size());
    for (int i = 0; i < residuals.size(); ++i)
        resvec(i) = residuals[i];
    
    octave_value_list retval;
    retval(0) = x;
    retval(1) = flag;
    retval(2) = relres;
    retval(3) = iter;
    retval(4) = resvec;
    return retval;
}
