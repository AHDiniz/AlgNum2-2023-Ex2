#include "octave/oct.h"
#include "octave/parse.h"
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
    octave_value_list retval(nargout);

    Matrix A = args(0).matrix_value();
    ColumnVector b = args(1).column_vector_value();
    double tol = args(2).double_value();
    int max_iter = args(3).int_value();
    
    const double EPS = octave::feval("eps")(0).double_value();

    if (tol >= 1.0)
        warning("Input tol is bigger than 1. \n Try to use a smaller tolerance.");
    else if (tol <= EPS / 2.0)
        warning("Input tol may not be achievable by pcg. \n Try to use a bigger tolerance.");

    int n = A.rows();
    double b_norm = euc_norm(b);

    ColumnVector x(n); // Solution
    for (int i = 0; i < n; ++i)
        x(i) = 0.0;
    ColumnVector x_pr(n); // Approximation in previous iteration
    ColumnVector x_min(n); // Approximation with minimum residual
    
    int flag = 0; // Convergence flag
    
    double relres = 0.0; // Relative residual by the end of the algorithm
    
    int iter = 0, iter_min = 0; // Number of iterations made and iteration with minimum residual
    
    ColumnVector resvec(max_iter); // Residuals' norms calculated in each iteration
    for (int i = 0; i < max_iter; ++i)
        resvec(i) = 0.0;

    if (b_norm == 0.0)
    {
        printf("The right hand side vector is all zero so pcg \n");
        printf("returned an all zero solution without iterating.\n");
        retval(0) = x;
        retval(1) = flag;
        retval(2) = relres;
        retval(3) = iter;
        retval(4) = resvec;
        return retval;
    }

    x_pr = x_min = x;

    ColumnVector r = b; // Residual vector
    resvec(0) = b_norm;
    
    // Auxiliar variables:
    double tau, old_tau, alpha, old_alpha;
    alpha = old_tau = 1.0;

    ColumnVector p(n);
    for (int i = 0; i < n; ++i)
        p(i) = 0.0;
    
    iter = 1;
    while (resvec(iter - 1) > tol * b_norm && iter < max_iter)
    {
        tau = dot_product(r, r);
        resvec(iter) = sqrt(tau);
        double beta = tau / old_tau;
        old_tau = tau;
        
        p = r + beta * p;
        
        ColumnVector w = A * p;
        
        old_alpha = alpha;
        alpha = tau / dot_product(p, w);

        x += alpha * p;
        r -= alpha * w;

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

    flag = (flag != 3 && iter >= max_iter) ? 1 : flag;

    retval(0) = x;
    retval(1) = flag;
    retval(2) = relres;
    retval(3) = iter;
    retval(4) = resvec;
    return retval;
}
