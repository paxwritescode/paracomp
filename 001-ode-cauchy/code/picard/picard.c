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

    

    /* Iterations of Picard method */
    do
    {
        for (j = 0; j < m; j++) 
        {
            double* yj = calloc(n, sizeof(double));
            double* yj1 = calloc(n, sizeof(double));

            for (i = 0; i < n; i++)
            {
                yj[i]  = y_prev[i][j];
                yj1[i] = y_prev[i][j + 1];
            }

            double *Ay_j  = matrix_mul_vector(n, n, A, yj);
            double *Ay_j1 = matrix_mul_vector(n, n, A, yj1);

            for (i = 0; i < n; i++)
            {
                y[i][j + 1] = y[i][j]
                    + delta_t / 2.0 * (
                        Ay_j[i]  + f[i][j]
                      + Ay_j1[i] + f[i][j + 1]
                    );
            }
            free(yj);
            free(yj1);
        }
        iterations++;
        for (i = 0; i < n; i++)
            for(j = 0; j < m; j++)
                y_prev[i][j] = y[i][j];
    } 
    while (compute_diff_norm(n, m + 1, y, y_prev) >= eps && iterations < 2);

    free_matrix(y_prev, n);
    printf("%d iterations\n", iterations);
    return y;
}