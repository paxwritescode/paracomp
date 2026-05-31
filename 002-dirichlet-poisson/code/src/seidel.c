#include "seidel.h"
#include "tools.h"
#include "def.h"

double **seidel(double border_x, double border_y,
                int Nx, int Ny,
                double eps,
                double (*f)(double, double),
                double (*u_left)(double), double (*u_right)(double),
                double (*u_up)(double), double (*u_down)(double),
                int max_iter)
{

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N_inner = Nx - 1;

    int base = N_inner / size;
    int rem  = N_inner % size;

    int local_Nx;
    int i_start, i_end;

    if (rank < rem) {
        local_Nx = base + 1;
        i_start = 1 + rank * local_Nx;
    } else {
        local_Nx = base;
        i_start = 1 + rem * (base + 1) + (rank - rem) * base;
    }

    i_end = i_start + local_Nx - 1;

    printf("rank %d: i_start=%d, i_end=%d, local_Nx=%d\n",
       rank, i_start, i_end, local_Nx);

    double **V_cur = alloc_matrix(local_Nx + 2, Ny + 1);
    double **V_prev = alloc_matrix(local_Nx + 2, Ny + 1);

    double h_x = border_x / (double)Nx;
    double h_y = border_y / (double)Ny;

    init_matrix(V_cur, local_Nx + 2, Ny + 1, 1.0);

    int i = 0, j = 0, iterations = 0;

    double diff = 1.0;

    /* LEFT boundary: i = 0 */
    if (rank == 0)
    {
        for (j = 0; j < Ny + 1; j++)
        {
            V_prev[0][j] = u_left(j * h_y);
            V_cur[0][j]  = V_prev[0][j];
        }
    }

    /* RIGHT boundary: i = Nx */
    if (rank == size - 1)
    {
        for (j = 0; j < Ny + 1; j++)
        {
            V_prev[local_Nx + 1][j] = u_right(j * h_y);
            V_cur[local_Nx + 1][j]  = V_prev[local_Nx + 1][j];
        }
    }

    for (i = 1; i <= local_Nx; i++)
    {
        int i_global = i_start + (i - 1);

        V_prev[i][0] = u_down(i_global * h_x);
        V_cur[i][0]  = V_prev[i][0];
    }

    for (i = 1; i <= local_Nx; i++)
    {
        int i_global = i_start + (i - 1);

        V_prev[i][Ny] = u_up(i_global * h_x);
        V_cur[i][Ny]  = V_prev[i][Ny];
    }

    double *send_left  = malloc((Ny + 1) * sizeof(double));
    double *send_right = malloc((Ny + 1) * sizeof(double));
    double *recv_left  = malloc((Ny + 1) * sizeof(double));
    double *recv_right = malloc((Ny + 1) * sizeof(double));

    /* CENTRAL NODES */
    do
    {
        int left  = rank - 1;
        int right = rank + 1;

        if (left < 0) left = MPI_PROC_NULL;
        if (right >= size) right = MPI_PROC_NULL;

        for (j = 0; j < Ny + 1; j++)
        {
            send_right[j] = V_prev[local_Nx][j];
            send_left[j]  = V_prev[1][j];
        }

        MPI_Sendrecv(
            send_right, Ny + 1, MPI_DOUBLE,
            right, 0,
            recv_left, Ny + 1, MPI_DOUBLE,
            left, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );

        MPI_Sendrecv(
            send_left, Ny + 1, MPI_DOUBLE,
            left, 1,
            recv_right, Ny + 1, MPI_DOUBLE,
            right, 1,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );

        if (left != MPI_PROC_NULL)
        {
            for (j = 0; j < Ny + 1; j++)
            {
                V_prev[0][j] = recv_left[j];
                V_cur[0][j]  = recv_left[j];
            }
        }
        if (right != MPI_PROC_NULL)
        {
            for (j = 0; j < Ny + 1; j++)
            {
                V_prev[local_Nx + 1][j] = recv_right[j];
                V_cur[local_Nx + 1][j]  = recv_right[j];
            }
        }

        for (i = 1; i <= local_Nx; i++)
        {
            int i_global = i_start + (i - 1);
            for (j = 1; j < Ny; j++)
                V_cur[i][j] =
                    ((((V_prev[i + 1][j] + V_cur[i - 1][j]) /
                       (h_x * h_x)) +
                      (V_prev[i][j + 1] + V_cur[i][j - 1]) /
                          (h_y * h_y)) -
                     f(i_global * h_x, j * h_y)) /
                    (2.0 / (h_x * h_x) + 2.0 / (h_y * h_y));
        }
        
        double local_diff = max_in_matrix_diff(V_prev, V_cur, local_Nx + 2, Ny + 1);
        MPI_Allreduce(
            &local_diff,
            &diff,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD
        );

        swap_matrices(&V_prev, &V_cur);

        iterations++;
    } while (diff > eps && iterations < max_iter);

    swap_matrices(&V_prev, &V_cur);

    free(send_left);
    free(send_right);
    free(recv_left);
    free(recv_right);

    free_matrix(V_prev, local_Nx + 2);

    if (rank == 0)
        printf("%d iterations\n", iterations);

    double *sendbuf = malloc(local_Nx * (Ny + 1) * sizeof(double));

    for (i = 1; i <= local_Nx; i++)
    {
        for (j = 0; j < Ny + 1; j++)
        {
            sendbuf[(i - 1)*(Ny + 1) + j] = V_cur[i][j];
        }
    }
    
    int local_size = local_Nx * (Ny + 1);

    int *recvcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));

    MPI_Gather(
        &local_size, 1, MPI_INT,
        recvcounts, 1, MPI_INT,
        0, MPI_COMM_WORLD
    );

    displs[0] = 0;
    for (int p = 1; p < size; p++)
        displs[p] = displs[p-1] + recvcounts[p-1];

    int total_size = displs[size-1] + recvcounts[size-1];

    double *recvbuf = NULL;
    if (rank == 0)
    {
        recvbuf = malloc(total_size * sizeof(double));
    }
    else
    {
        recvbuf = malloc(1 * sizeof(double));
    }

    MPI_Gatherv(
        sendbuf, local_size, MPI_DOUBLE,
        recvbuf, recvcounts, displs, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    double **V_global = NULL;

    if (rank == 0)
    {
        V_global = alloc_matrix(Nx + 1, Ny + 1);

        int offset = 0;

        for (int p = 0; p < size; p++)
        {
            int p_local_Nx = (p < rem) ? base + 1 : base;

            int p_i_start = (p < rem)
                ? 1 + p * (base + 1)
                : 1 + rem * (base + 1) + (p - rem) * base;

            for (int i_loc = 0; i_loc < p_local_Nx; i_loc++)
            {
                int i_global = p_i_start + i_loc;

                for (j = 0; j < Ny + 1; j++)
                {
                    V_global[i_global][j] =
                        recvbuf[offset + i_loc*(Ny+1) + j];
                }
            }

            offset += recvcounts[p];
        }


        for (j = 0; j < Ny + 1; j++)
        {
            V_global[0][j] = u_left(j * h_y);
            V_global[Nx][j] = u_right(j * h_y);
        }

        if (Nx == 80 && Ny == 400)
        {
            FILE *f_sol = fopen("results/solution.csv", "w");
            if (f_sol != NULL)
            {
                for (int gi = 0; gi <= Nx; gi++)
                {
                    for (int gj = 0; gj <= Ny; gj++)
                    {
                        fprintf(f_sol, "%lf%s", V_global[gi][gj], (gj == Ny) ? "" : ",");
                    }
                    fprintf(f_sol, "\n");
                }
                fclose(f_sol);
                printf("[Rank 0] Numerical solution successfully saved to 'results/solution.csv'\n");
            }
            else
            {
                printf("[Rank 0] Error opening results/solution.csv for writing!\n");
            }
        }

        free_matrix(V_cur, local_Nx + 2);
        free(sendbuf);
        free(recvbuf);
        free(recvcounts);
        free(displs);

        return V_global;
    }

    free_matrix(V_cur, local_Nx + 2);
    free(sendbuf);
    free(recvbuf);
    free(recvcounts);
    free(displs);

    return NULL;
}