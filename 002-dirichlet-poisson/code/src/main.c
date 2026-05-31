#include "def.h"
#include "test.h"
#include "study.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // run_test_case(rank);
    // study_time_vs_grid_size(rank, size);

    // study_scaling_performance(rank, size);

    double **V = seidel(PI, 1.0, 
                        80, 400, 1e-6, 
                        &f, &u_0y, &u_piy, &u_x1, &u_x0, 
                        100000);

    if (V) free_matrix(V, 80 + 1);

    MPI_Finalize();

    return 0;
}