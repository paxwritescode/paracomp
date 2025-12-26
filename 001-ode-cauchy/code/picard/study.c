#include "study.h"

void study_time_vs_size(void)
{
    double t0 = 0.0;
    double t = 1.0;
    double eps = 1e-6;

    int k = 0;
    int cases = 20;

    int n_array[cases];
    int m = 1000;
    for (k = 0; k < cases; k++)
    {
        n_array[k] = 50 * (k + 1);
        printf("n_array[%d] = %d\n", k, n_array[k]);
    }

    // int p = omp_get_max_threads();
    int p = 4;
    omp_set_num_threads(p);

    FILE *f_file = fopen("results/time_vs_size.csv", "w");
    if (!f_file)
    {
        perror("fopen");
        return;
    }

    printf("Time vs size\n");
    printf("%d threads\n\n", p);

    fprintf(f_file, "# n = 10 * k, k from 1 to 25\n\n");

    for (k = 0; k < cases; k++)
    {
        int n = n_array[k];

        printf("k = %d, n_array[%d] = %d\n", k, k, n);
        double **A = generate_matrix(n, 1.0 / n);

        double **f = generate_rhs(n, m, t0, t);

        double t_start = omp_get_wtime();
        double **y = picard_method(n, f, A, eps, t0, t, m);
        double t_end = omp_get_wtime();

        printf("%d, %lf\n", n, t_end - t_start);
        fprintf(f_file, "%.6lf, ", t_end - t_start);

        free_matrix(y, n);
        free_matrix(f, n);
        free_matrix(A, n);
    }

    fclose(f_file);
}

void study_dependencies_on_threads(void)
{
    double t0 = 0.0;
    double t = 1.0;
    double eps = 1e-6;

    int n = 300;
    int m = 3000;

    int k = 0;

    int max_threads = omp_get_max_threads();
    int threads[32];
    int num_cases = generate_threads(threads, max_threads);

    FILE *file_time = fopen("results/time_vs_threads.csv", "w");
    FILE *file_speedup = fopen("results/speedup_vs_threads.csv", "w");
    FILE *file_efficiency = fopen("results/efficiency_vs_threads.csv", "w");

    if (!file_time || !file_speedup || !file_efficiency)
    {
        perror("fopen");
        return;
    }

    fprintf(file_time, "# threads from 1 to 12\n\n");
    fprintf(file_speedup, "# threads from 1 to 12\n\n");
    fprintf(file_efficiency, "# threads from 1 to 12\n\n");

    double **A = generate_matrix(n, 1.0 / n);
    double **rhs = generate_rhs(n, m, t0, t);

    /* 1 thread, baseline run */
    omp_set_num_threads(1);
    printf("Baseline run (1 thread)\n");
    double t_start = omp_get_wtime();
    double **y = picard_method(n, rhs, A, eps, t0, t, m);
    double t_end = omp_get_wtime();

    double T1 = t_end - t_start;

    fprintf(file_time, "%.6lf, ", T1);
    fprintf(file_speedup, "1.000000, ");
    fprintf(file_efficiency, "1.000000, ");
    free_matrix(y, n);

    /* Parallel runs */
    for (k = 0; k < num_cases; k++)
    {
        int p = threads[k];
        if (p == 1)
            continue;

        omp_set_num_threads(p);
        printf("p = %d\n", p);

        t_start = omp_get_wtime();
        y = picard_method(n, rhs, A, eps, t0, t, m);
        t_end = omp_get_wtime();

        double Tp = t_end - t_start;
        double speedup = T1 / Tp;
        double efficiency = T1 / (p * Tp);

        fprintf(file_time, "%.6lf, ", Tp);
        fprintf(file_speedup, "%.6lf, ", speedup);
        fprintf(file_efficiency, "%.6lf, ", efficiency);

        free_matrix(y, n);
    }

    free_matrix(rhs, n);
    free_matrix(A, n);

    fclose(file_time);
    fclose(file_speedup);
    fclose(file_efficiency);
}
