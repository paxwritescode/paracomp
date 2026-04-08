#include "test.h"

void run_test_case(int rank)
{
    if (rank == 0)
        printf("TEST CASE FROM REPORT: \n\n");

    int Nx_test = 4, Ny_test = 2, iter_test = 2;
    double right_border = PI, upper_border = 1.0;

    double eps = 1e-6;

    double **V = seidel(right_border, upper_border, Nx_test, Ny_test, eps, &f, &u_0y, &u_piy, &u_x1, &u_x0, iter_test);

    if (rank == 0)
    {
        for (int j = Ny_test; j >= 0; j--)
        {
            for (int i = 0; i < Nx_test + 1; i++)
                printf("%.2lf ", V[i][j]);
            printf("\n");
        }
        free_matrix(V, Nx_test + 1);
    }
}