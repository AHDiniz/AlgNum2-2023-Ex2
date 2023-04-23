#include "octave/oct.h"
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

    int n = A.rows();

    ColumnVector x0(n);
    for (int i = 0; i < n; ++i)
        x0(i) = 0.0;
    ColumnVector d0 = b;
    ColumnVector r0 = b;
    ColumnVector x(n);
    ColumnVector d(n);
    ColumnVector r(n);

    double delta = dot_product(r0, r0);
    double delta0 = delta;
    double epsilon = tol * tol * delta0;

    int iter = 0;
    int flag = 0;

    ColumnVector resvec(max_iter);
    for (int i = 0; i < max_iter; ++i)
        resvec(i) = 0.0;

    while (delta < epsilon && iter < max_iter)
    {
        ColumnVector v = A * d0;
        double lambda = delta / dot_product(d0, v);
        x = x0 + lambda * d0;
        r = r0 - lambda * v;
        resvec(iter) = euc_norm(r);
        delta0 = delta;
        delta = dot_product(r, r);
        double beta = delta0 / delta;
        d = r + beta * d0;
        x0 = x;
        r0 = r;
        d0 = d;
        ++iter;

        if (euc_norm(x - x0) <= epsilon * euc_norm(x))
        {
            flag = 3;
            break;
        }
    }

    if (iter >= max_iter)
        flag = 1;

    octave_value_list retval(nargout);
    retval(0) = x;
    retval(1) = flag;
    retval(2) = resvec(iter - 1);
    retval(3) = iter;
    retval(4) = resvec;
    return retval;
}
