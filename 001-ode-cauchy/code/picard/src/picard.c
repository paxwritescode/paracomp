#include <omp.h>

#include "picard.h"
#include "generate.h"
#include "matrix_tools.h"

double **picard_method(int n, double **f, double **A, double eps, double t_0, double t, int m, FILE *f_diff, FILE *f_norm)
{
    int i = 0, j = 0;
    int iterations = 0;

    double delta_t = (t - t_0) / (double)m;

    double **y_prev = alloc_matrix(n, m + 1);
    double **y = alloc_matrix(n, m + 1);

    double diff = 0;

    double initial = 1.0;

    double solution_norm_j = 0;

    /* Initial condition */
    for (i = 0; i < n; i++)
        for (j = 0; j < m + 1; j++)
            y_prev[i][j] = 1.0;

    /* Iterations of Picard method */
    do
    {
        for (i = 0; i < n; i++)
            y[i][0] = 1.0;

#pragma omp parallel private(i)
        {
            double *yj = calloc(n, sizeof(double));
            double *yj1 = calloc(n, sizeof(double));

            double *Ay_j = calloc(n, sizeof(double));
            double *Ay_j1 = calloc(n, sizeof(double));

#pragma omp for schedule(static)
            for (j = 0; j < m; j++)
            {
                for (i = 0; i < n; i++)
                {
                    yj[i] = y_prev[i][j];
                    yj1[i] = y_prev[i][j + 1];
                }

                matrix_mul_vector(n, A, yj, Ay_j);
                matrix_mul_vector(n, A, yj1, Ay_j1);

                for (i = 0; i < n; i++)
                {
                    y[i][j + 1] = initial + delta_t / 2.0 * (Ay_j[i] + f[i][j] + Ay_j1[i] + f[i][j + 1]);
                }
            }
            free(yj);
            free(yj1);
            free(Ay_j);
            free(Ay_j1);
        }

        diff = compute_diff_norm(n, m + 1, y, y_prev);
        if (f_diff)
            fprintf(f_diff, "%.6lf, ", diff);

        swap_arrays(&y, &y_prev);

        iterations++;

        printf("Iter: diff = %.6e\n", diff);

    } while (diff > eps);


    printf("%d iterations\n", iterations);

    for (j = 0; j <= m; j++)
    {
        double sum = 0.0;
        for (i = 0; i < n; i++)
        {
            sum += y[i][j] * y[i][j];
        }

        solution_norm_j = sqrt(sum);
        if (f_norm)
            fprintf(f_norm, "%.6lf, ", solution_norm_j);
    }

    free_matrix(y_prev, n);
    return y;
}