#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "picard.h"
#include "generate.h"
#include "matrices.h"

double **picard_method(int n, double **f, double**A, double eps, double t_0, double t, int m)
{
    int i = 0, j = 0;
    int iterations = 0;

    double delta_t = (t - t_0) / (double)m;

    double **y_prev = alloc_matrix(n, m + 1);
    double **y = alloc_matrix(n, m + 1);

    double diff = 0;

    double* yj = calloc(n, sizeof(double));
    double* yj1 = calloc(n, sizeof(double));

    double* Ay_j = calloc(n, sizeof(double));
    double* Ay_j1 = calloc(n, sizeof(double));

    /* Iterations of Picard method */
    do
    {
        for (j = 0; j < m; j++) 
        {
            for (i = 0; i < n; i++)
            {
                yj[i]  = y_prev[i][j];
                yj1[i] = y_prev[i][j + 1];
            }

            matrix_mul_vector(n, A, yj, Ay_j);
            matrix_mul_vector(n, A, yj1, Ay_j1);

            for (i = 0; i < n; i++)
            {
                y[i][j + 1] = y_prev[i][j]
                    + delta_t / 2.0 * (
                        Ay_j[i]  + f[i][j]
                      + Ay_j1[i] + f[i][j + 1]
                    );
            }
        }

        diff = compute_diff_norm(n, m + 1, y, y_prev);

        for (i = 0; i < n; i++)
            for(j = 0; j < m + 1; j++)
                y_prev[i][j] = y[i][j];
        iterations++;
    } 
    while (diff > eps);

    free(yj); free(yj1);
    free(Ay_j); free(Ay_j1);

    free_matrix(y_prev, n);

    printf("%d iterations\n", iterations);

    return y;
}