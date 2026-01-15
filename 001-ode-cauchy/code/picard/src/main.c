#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "matrix_tools.h"
#include "generate.h"
#include "picard.h"
#include "study.h"

void run_test_case(void)
{
    int n_0 = 3;
    int m_0 = 3;
    double t_start = 0.0;
    double t_end = 0.3;
    double eps = 1e-5;

    /* test matrix A  */
    printf("TEST CASE FROM REPORT\n\n");
    double **A = alloc_matrix(n_0, n_0);
    for (int i = 0; i < n_0; i++)
        for (int j = 0; j < n_0; j++)
            A[i][j] = (i <= j) ? 1.0 : 0.0;

    /* test RHS */
    // double **f = generate_test_rhs(n, m, t_0, t);
    double **f = generate_rhs(n_0, m_0, t_start, t_end);

    /* run Picard */
    double **y = picard_method(n_0, f, A, eps, t_start, t_end, m_0, NULL, NULL);

    /* print result */
    printf("Result y:\n");
    for (int i = 0; i < n_0; i++)
    {
        for (int j = 0; j <= m_0; j++)
            printf("%8.5f ", y[i][j]);
        printf("\n");
    }

    /* cleanup */
    free_matrix(y, n_0);
    free_matrix(f, n_0);
    free_matrix(A, n_0);
}

int main(void)
{
    // run_test_case();

    // study_time_vs_size();
    study_dependencies_on_threads(250.0);
    print_norm_diff(3.0);

    return 0;
}
