#include "seidel.h"
#include "tools.h"

void exchange_boundaries(double **V, int local_Nx, int Ny, int rank, int size)
{
    MPI_Status status;

    if (rank < size - 1)
    {
        MPI_Sendrecv(
            V[local_Nx], Ny + 1, MPI_DOUBLE, rank + 1, 0,
            V[local_Nx + 1], Ny + 1, MPI_DOUBLE, rank + 1, 0,
            MPI_COMM_WORLD, &status
        );
    }

    if (rank > 0)
    {
        MPI_Sendrecv(
            V[1], Ny + 1, MPI_DOUBLE, rank - 1, 0,
            V[0], Ny + 1, MPI_DOUBLE, rank - 1, 0,
            MPI_COMM_WORLD, &status
        );
    }
}

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

    // local mesh initialization
    for (int i_local = 0; i_local <= local_Nx + 1; i_local++)
    {
        int i_global = start_x + i_local - 1;

        if (i_global >= 0 && i_global <= Nx)
        {
            for (int j = 0; j <= Ny; j++)
            {
                // lower
                if (j == 0)
                    V[i_local][j] = u_down(i_global * h_x);
                // upper
                else if (j == Ny)
                    V[i_local][j] = u_up(i_global * h_x);
                else
                    V[i_local][j] = 0.0; // start value
            }
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

    // main cycle
    do
    {
        diff = 0.0;

        // RED UPDATE

        for (int i = 1; i <= local_Nx; i++) 
        {
            int i_global = start_x + i - 1;
            if (i_global < 0 || i_global > Nx)
            {
                printf("Error in i_global\n");
                return NULL;
            }

            if (i_global == 0 || i_global == Nx)
                continue;

            int j_start = (i_global % 2 == 0) ? 1 : 2;

            for (int j = j_start; j < Ny; j += 2)
            {
                double old = V[i][j];

                V[i][j] =
                    (((V[i + 1][j] + V[i - 1][j]) / (h_x * h_x)) +
                    ((V[i][j + 1] + V[i][j - 1]) / (h_y * h_y)) -
                    f(i_global * h_x, j * h_y)) /
                    (2.0 / (h_x * h_x) + 2.0 / (h_y * h_y));
            
                double local_diff = fabs(V[i][j] - old);
                if (local_diff > diff)
                    diff = local_diff;
            }
        }

        exchange_boundaries(V, local_Nx, Ny, rank, size);

        // BLACK UPDATE

        for (int i = 1; i <= local_Nx; i++)
        {
            int i_global = start_x + i - 1;

            if (i_global == 0 || i_global == Nx)
                continue;

            int j_start = (i_global % 2 == 0) ? 2 : 1;

            for (int j = j_start; j < Ny; j += 2)
            {
                double old = V[i][j];

                V[i][j] =
                    (((V[i + 1][j] + V[i - 1][j]) / (h_x * h_x)) +
                    ((V[i][j + 1] + V[i][j - 1]) / (h_y * h_y)) -
                    f(i_global * h_x, j * h_y)) /
                    (2.0 / (h_x * h_x) + 2.0 / (h_y * h_y));

                double local_diff = fabs(V[i][j] - old);
                if (local_diff > diff)
                    diff = local_diff;
            }
        }

        exchange_boundaries(V, local_Nx, Ny, rank, size);
    
        double global_diff;

        MPI_Allreduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        diff = global_diff;

        iterations++;
        // if (rank == 0)
        //     printf("iter %d global diff = %.10lf\n", iterations, diff);
    } while (diff > eps && iterations < max_iter);
    
    int send_count = local_Nx * (Ny + 1);

    double *sendbuf = malloc(send_count * sizeof(double));

    int k = 0;
    for (int i = 1; i <= local_Nx; i++)
        for (int j = 0; j <= Ny; j++)
            sendbuf[k++] = V[i][j];

    double *recvbuf = NULL;
    int *recvcounts = NULL;
    int *displs = NULL;

    if (rank == 0)
    {
        recvbuf = malloc((Nx + 1) * (Ny + 1) * sizeof(double));
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
    }

    int local_count = send_count;

    if (rank == 0)
    {
        int offset = 0;
        for (int p = 0; p < size; p++)
        {
            int p_local_Nx = (p < remainder) ? base + 1 : base;

            recvcounts[p] = p_local_Nx * (Ny + 1);
            displs[p] = offset;

            offset += recvcounts[p];
        }
    }

    MPI_Gatherv(
        sendbuf,
        send_count,
        MPI_DOUBLE,
        recvbuf,
        recvcounts,
        displs,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );
    free(sendbuf);

    double **V_global = NULL;

    if (rank == 0)
    {
        V_global = alloc_matrix(Nx + 1, Ny + 1);

        int k = 0;
        for (int i = 0; i <= Nx; i++)
            for (int j = 0; j <= Ny; j++)
                V_global[i][j] = recvbuf[k++];
    }

    free_matrix(V, local_width);

    if (rank == 0)
    {
        free(recvbuf);
        free(recvcounts);
        free(displs);
    }

    if (rank == 0)
    {
        printf("%d iterations\n", iterations);
        return V_global;
    }
    else
        return NULL;
}