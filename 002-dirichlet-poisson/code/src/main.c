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

    study_scaling_performance(rank, size);

    MPI_Finalize();

    return 0;
}