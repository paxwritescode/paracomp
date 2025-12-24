#include "study.h"

void study_time_vs_size(void)
{
    double t0 = 0.0;
    double t = 1.0;
    double eps = 1e-6;

    int k = 0;
    int cases = 25;

    int n_array[cases];
    int m = 100;
    for (k = 1; k <= cases; k++)
        n_array[k] = 10 * k;

    int p = omp_get_max_threads();
    omp_set_num_threads(p);

    FILE* f_file = fopen("results/time_vs_size.csv", "w");

    printf("%d threads\n\n", p);
    printf("n, t\n");

    fprintf(f_file, "# n = 10 * k, k from 1 to 25\n\n");

    for (k = 0; k < cases; k++)
    {
        int n = n_array[k];

        double **A = alloc_matrix(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = (i >= j) ? 1.0 : 0.0;

        double** f = generate_rhs(n, m, t0, t);

        double t_start = omp_get_wtime();
        double **y = picard_method(n, f, A, eps, t0, t, m);
        double t_end = omp_get_wtime();

        //printf("%d, %lf\n", n, t_end - t_start);
        fprintf(f_file, "%.6lf, ", t_end - t_start);

        free_matrix(y, n);
        free_matrix(f, n);
        free_matrix(A, n);
    }

    fclose(f_file);
}