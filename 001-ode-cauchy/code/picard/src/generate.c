#include "generate.h"
#include "matrix_tools.h"

/* n dimensions, m + 1 nodes */

double **generate_rhs(int n, int m, double t0, double t)
{
    double **f = alloc_matrix(n, m + 1);

    int i, j = 0;
    double delta_t = (t - t0) / (double)(m);
    for (i = 0; i < n; i++)
        for (j = 0; j < m + 1; j++)
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
        p += 1;
    }
    return k;
}

double **generate_matrix(int n)
{
    double **A = alloc_matrix(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = (i >= j) ? 1.0 : 0.0;

    return A;
}