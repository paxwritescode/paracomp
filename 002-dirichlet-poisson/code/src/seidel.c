#include "seidel.h"
#include "tools.h"

double **seidel(double border_x, double border_y,
                int Nx, int Ny,
                double eps,
                double (*f)(double, double),
                double (*u_left)(double), double (*u_right)(double),
                double (*u_up)(double), double (*u_down)(double),
                int max_iter)
{
    double **V_cur = alloc_matrix(Nx + 1, Ny + 1);
    double **V_prev = alloc_matrix(Nx + 1, Ny + 1);

    double h_x = border_x / (double)Nx;
    double h_y = border_y / (double)Ny;

    init_matrix(V_cur, Nx + 1, Ny + 1, 1.0);

    int i = 0, j = 0, iterations = 0;

    double diff = 1.0;

    /* BORDER NODES */

    /* left & right */
    for (j = 0; j < Ny + 1; j++)
    {
        V_prev[0][j] = u_left(j * h_y);
        V_cur[0][j] = V_prev[0][j];

        V_prev[Nx][j] = u_right(j * h_y);
        V_cur[Nx][j] = V_prev[Nx][j];
    }

    /* lower & upper */
    for (i = 0; i < Nx + 1; i++)
    {
        V_prev[i][0] = u_down(i * h_x);
        V_cur[i][0] = V_prev[i][0];

        V_prev[i][Ny] = u_up(i * h_x);
        V_cur[i][Ny] = V_prev[i][Ny];
    }

    /* CENTRAL NODES */
    do
    {
        for (i = 1; i < Nx; i++)
            for (j = 1; j < Ny; j++)
                V_cur[i][j] =
                    ((((V_prev[i + 1][j] + V_cur[i - 1][j]) /
                       (h_x * h_x)) +
                      (V_prev[i][j + 1] + V_cur[i][j - 1]) /
                          (h_y * h_y)) -
                     f(i * h_x, j * h_y)) /
                    (2.0 / (h_x * h_x) + 2.0 / (h_y * h_y));
        
        diff = max_in_matrix_diff(V_prev, V_cur, Nx + 1, Ny + 1);

        swap_matrices(&V_prev, &V_cur);

        iterations++;
    } while (diff > eps && iterations < max_iter);

    swap_matrices(&V_prev, &V_cur);

    free_matrix(V_prev, Nx + 1);
    printf("%d iterations\n", iterations);

    return V_cur;
}