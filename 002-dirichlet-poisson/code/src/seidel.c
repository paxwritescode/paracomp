#include "seidel.h"
#include "tools.h"

double **seidel(double border_x, double border_y,
                int Nx, int Ny,
                double eps,
                double (*f)(double, double),
                double (*u_left)(double), double (*u_right)(double),
                double (*u_up)(double), double (*u_down)(double),
                int max_iter,
                int size, int rank)
{
    if (size > Nx)
    {
        printf("Domain size is not valid: size = %d, Nx = %d\n", size, Nx);
        return NULL;
    }

    int local_Nx;

    double h_x = border_x / (double)Nx;
    double h_y = border_y / (double)Ny;

    int i = 0, j = 0, iterations = 0;

    double diff = 1.0;

    // (Nx + 1) nodes
    int base = (Nx + 1) / size;
    int remainder = (Nx + 1) % size;

    // distribution
    if (rank < remainder)
        local_Nx = base + 1;
    else
        local_Nx = base;

    // + 2 ghost columns
    int local_width = local_Nx + 2;
    double **V = alloc_matrix(local_width, Ny + 1);

    int start_x = rank * base + (rank < remainder ? rank : remainder);

    // diagnostic print
    for (int p = 0; p < size; p++)
    {
        if (rank == p)
        {
            printf("rank %d: local_Nx = %d, start_x = %d\n\n",
                rank, local_Nx, start_x);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // local mesh initialization
    for (int i_local = 0; i_local < local_Nx + 2; i_local++)
    {
        int i_global = start_x + i_local - 1; 

        for (int j = 0; j <= Ny; j++)
        {
            // lower
            if (j == 0)
                V[i_local][j] = u_down(i_global * h_x);
            // upper
            else if (j == Ny)
                V[i_local][j] = u_up(i_global * h_x);
            else
                V[i_local][j] = 1.0; // start value
        }
    }

    if (rank == 0)
    {
        for (int j = 0; j <= Ny; j++)
            V[1][j] = u_left(j * h_y);
    }

    if (rank == size - 1)
    {
        for (int j = 0; j <= Ny; j++)
            V[local_Nx][j] = u_right(j * h_y); 
    }

    for (int p = 0; p < size; p++)
    {
        if (rank == p)
        {
            printf("rank %d local grid:\n", rank);
            for (int j = Ny; j >= 0; j--)
            {
                for (int i = 0; i <= local_Nx + 1; i++)
                    printf("%.4lf ", V[i][j]);
                printf("\n");
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    return NULL;
}