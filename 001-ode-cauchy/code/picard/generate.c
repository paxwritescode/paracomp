#include "generate.h"
#include "matrix_tools.h"

/* n dimensions, m + 1 nodes */

/* test case: t_0 = 0, t = 0.3, m = 3, n = 3*/
double **generate_test_rhs(int n, int m, double t_0, double t)
{
    double delta_t = (t - t_0) / m;
    double **f = alloc_matrix(n, m + 1);

    for (int j = 0; j <= m; j++)
    {
        double tj = t_0 + j * delta_t;

        f[0][j] = tj;
        f[1][j] = tj * tj;
        f[2][j] = sin(tj);
    }

    return f;
}


double** generate_rhs(int n, int m, double t0, double t)
{
    double** f = alloc_matrix(n, m + 1);

    int i, j = 0;
    double delta_t = (t - t0) / (double)(m);
    for (i = 0; i < n; i++)
        for(j = 0; j < m + 1; j++)
            f[i][j] = sin(t0 + j * delta_t + i);

    return f;
}

int generate_threads(int *threads, int max_threads)
{
    int k = 0;
    int p = 1;

    while (p <= max_threads)
    {
        threads[k++] = p;
        p *= 2;
    }
    return k;
}