#include "seidel.h"
#include "tools.h"

double** seidel(double border_x, double border_y, int Nx, int Ny, double eps)
{
    double** V_cur = alloc_matrix((Nx - 1), (Ny - 1));
    double** V_prev = alloc_matrix((Nx - 1), (Ny - 1));

    double h_x = border_x / (double)Nx;
    double h_y = border_y / (double)Ny;

    while (max_in_matrix_diff(V_prev, V_cur, (Nx - 1), (Ny - 1)) > eps)
    {
        /* code */
    }
    
    free(V_prev);
    return V_cur;
}