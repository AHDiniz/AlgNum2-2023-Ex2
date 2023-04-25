#include "octave/oct.h"
#include "octave/parse.h"
#include <vector>
#include <cmath>

inline double dot_product(const ColumnVector &a, const ColumnVector &b)
{
    double dot = 0.0;
    for (int i = 0; i < a.numel(); ++i)
        dot += a(i) * b(i);
    return dot;
}

inline double euc_norm(const ColumnVector &a)
{
    return sqrt(dot_product(a, a));
}

inline ColumnVector gauss(const Matrix &M, const ColumnVector &v)
{
    octave_value_list gaussIN;
    gaussIN(0) = M;
    gaussIN(1) = v;
    return octave::feval("\\", gaussIN)(0).column_vector_value();
}

DEFUN_DLD(pc_conjgrad, args, nargout, "C++ implementation of the Conjugate Gradient method.")
{
    Matrix A, M1, M2;
    ColumnVector b;
    double tol;
    int maxit;
    bool precond = false;

    A = args(0).matrix_value();
    b = args(1).column_vector_value();
    tol = args(2).double_value();
    maxit = args(3).int_value();

    int n = A.rows();
    ColumnVector x(n);
    for (int i = 0; i < n; ++i)
        x(i) = 0;
    int flag = 0;
    double relres = 0.0;
    int iter = 0;
    ColumnVector resvec(maxit + 2);
    for (int i = 0; i < maxit + 2; ++i)
        resvec(i) = 0;

    if (args.length() > 4)
    {
        M1 = args(4).matrix_value();
        M2 = args(5).matrix_value();
        precond = true;
    }

    const double EPS = octave::feval("eps")(0).double_value();

    if (tol >= 1)
        warning("Input tol is bigger then 1.\nTry using a smaller tolerance.");
    else if (tol <= EPS / 2)
        warning("Input tol may not be achievable by pc_conjgrad.\nTry using a bigger tolerance.");
    
    maxit += 2;

    double bNorm = euc_norm(b);
    if (bNorm == 0.0)
    {
        octave_value_list retval;
        retval(0) = x;
        retval(1) = flag;
        retval(2) = relres;
        retval(3) = iter;
        retval(4) = resvec;
        return retval;
    }

    ColumnVector x_pr = x;
    ColumnVector x_min = x;

    ColumnVector r = b;
    iter = 1;
    int iter_min = 0;
    flag = 1;
    resvec(0) = bNorm;
    ColumnVector p = x;
    double alpha = 1.0, old_tau = 1.0;

    while (resvec(iter - 1) > tol * bNorm && iter < maxit)
    {
        ColumnVector z = r;
        if (precond)
        {
            z = gauss(M1, r);
            z = gauss(M2, z);
        }

        double tau = dot_product(z, r);
        resvec(iter - 1) = sqrt(tau);
        double beta = tau / old_tau;
        old_tau = tau;
        p = z + beta * p;
        ColumnVector w = A * p;

        double den = dot_product(p, w);
        alpha = tau / den;

        x += alpha * p;
        r -= alpha * w;
        resvec(iter) = euc_norm(r);

        if (resvec(iter) <= resvec(iter_min))
        {
            x_min = x;
            iter_min = iter;
        }

        ++iter;
        if (euc_norm(x - x_pr) <= EPS * euc_norm(x))
        {
            flag = 3;
            break;
        }

        x_pr = x;
    }

    ColumnVector true_resvec(iter - 1);
    for (int i = 0; i < iter - 1; ++i)
        true_resvec(i) = resvec(i);
    resvec = true_resvec;

    relres = (resvec(0) == 0.0) ? 0.0 : resvec(iter_min + 1) / bNorm;

    iter -= 2;

    if (flag == 1 && relres <= tol)
        flag = 0;

    octave_value_list retval;
    retval(0) = x;
    retval(1) = flag;
    retval(2) = relres;
    retval(3) = iter;
    retval(4) = resvec;
    return retval;
}