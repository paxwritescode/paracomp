#include <stdio.h>
#include <mpi.h>

#include "seidel.h"
#include "def.h"
#include "functions.h"
#include "tools.h"

void run_test_case(int rank, int size)
{
    if (rank == 0)
        printf("TEST CASE FROM REPORT: \n\n");

    int Nx_test = 4, Ny_test = 2, iter_test = 2;
    double right_border = PI, upper_border = 1.0;

    double eps = 1e-6;

    int local_Nx = 0;

    double **V = seidel(right_border, upper_border, 
        Nx_test, Ny_test, 
        eps, 
        &f, 
        &u_0y, &u_piy, &u_x1, &u_x0, 
        iter_test, 
        size, rank);

    if (V != NULL)
    {
        printf("\nSolution:\n");
        for (int j = Ny_test; j >= 0; j--)
        {
            for (int i = 0; i < Nx_test + 1; i++)
            {
                printf("%.4lf ", V[i][j]);
            }
            printf("\n");
        }

        free_matrix(V, Nx_test + 1);
    }

}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    run_test_case(rank, size);

    MPI_Finalize();

    return 0;
}