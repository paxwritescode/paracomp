#include "test.h"

void run_test_case(int rank)
{
    if (rank == 0)
        printf("TEST CASE FROM REPORT: \n\n");

    int Nx_test = 200, Ny_test = 200, iter_test = 200000;
    double right_border = PI, upper_border = 1.0;

    double eps = 1e-6;

    double t_start = MPI_Wtime();
    double **V = seidel(right_border, upper_border, 
        Nx_test, Ny_test, 
        eps, 
        &f, 
        &u_0y, &u_piy, &u_x1, &u_x0, 
        iter_test);

    double t_end = MPI_Wtime();
    double local_time = t_end - t_start;

    double global_time;

    MPI_Reduce(&local_time, &global_time, 1,
            MPI_DOUBLE, MPI_MAX,
            0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Time: %lf seconds\n", global_time);
        free_matrix(V, Nx_test + 1);
    }

}
