#include "study.h"

#include <stdio.h>
#include <mpi.h>
#include "seidel.h"
#include "functions.h"
#include "tools.h"

void study_dependencies_on_nodes(void)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 4)
    {
        if (rank == 0)
            printf("This study assumes exactly 4 MPI processes.\n");
        return;
    }

    FILE *f_time_on_nodes = fopen("results/time_on_nodes.csv", "w");
    if (f_time_on_nodes == NULL)
    {
        printf("Error: file not opened\n");
        return;
    }

    int Nxy_max = 256;
    fprintf(f_time_on_nodes, "# Nx and Ny simultaneously from 4 to %d, size = 4\n", Nxy_max);
    fprintf(f_time_on_nodes, "# 2run: mpirun -np 4 -env FI_PROVIDER=tcp ./build/seidel \n\n");


    double border_x = PI;
    double border_y = 1.0;
    double eps = 1e-5;
    int max_iter = 10000000;


    int Ny = 4;
    for (int Nx = 4; Nx <= Nxy_max; Nx *= 2)
    {
        MPI_Barrier(MPI_COMM_WORLD); 
        double t_start = MPI_Wtime();

        double **V_global = seidel(border_x, border_y,
                                   Nx, Ny,
                                   eps,
                                   &f,
                                   &u_0y, &u_piy, &u_x1, &u_x0, 
                                   max_iter,
                                   size, rank);

        MPI_Barrier(MPI_COMM_WORLD);
        double t_end = MPI_Wtime();

        if (rank == 0)
        {
            printf("Nx = %d, Ny = %d -> Time: %.6lf s\n", Nx, Ny, t_end - t_start);
            fprintf(f_time_on_nodes, "%.6lf, ", t_end - t_start);
        }

        if (rank == 0 && V_global != NULL)
            free_matrix(V_global, Nx + 1);

        Ny *= 2;
    }

    fclose(f_time_on_nodes);
}

void study_dependencies_on_threads(int rank, int size)
{
    int Nx = 64;
    int Ny = 64;
    int max_iter = 10000;
    double eps = 1e-6;

    double border_x = PI;
    double border_y = 1.0;

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    seidel(border_x, border_y,
           Nx, Ny,
           eps,
           &f,
           &u_0y, &u_piy, &u_x1, &u_x0, 
           max_iter,
           size,
           rank);

    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

    if (rank == 0)
        printf("Time = %lf\n", end - start);
}