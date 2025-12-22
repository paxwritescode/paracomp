#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "matrix_tools.h"
#include "generate.h"
#include "picard.h"

void run_test_case(void)
{
    int n = 3;
    int m = 3;
    double t_0 = 0.0;
    double t = 0.3;
    double eps = 1e-6;

    /* test matrix A  */
    double **A = alloc_matrix(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = (i >= j) ? 1.0 : 0.0;

    /* test RHS */
    double **f = generate_test_rhs(n, m, t_0, t);

    printf("TEST CASE FROM REPORT\n\n");

    /* run Picard */
    double **y = picard_method(n, f, A, eps, t_0, t, m);

    /* print result */
    printf("Result y:\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= m; j++)
            printf("%8.5f ", y[i][j]);
        printf("\n");
    }

    /* cleanup */
    free_matrix(y, n);
    free_matrix(f, n);
    free_matrix(A, n);
}

int main(void)
{
    run_test_case();
    // #pragma omp parallel
    // {
    //     printf("thread %d\n", omp_get_thread_num());
    // }
    return 0;
}
