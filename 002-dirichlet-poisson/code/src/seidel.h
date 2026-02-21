#ifndef SEIDEL_H
#define SEIDEL_H

#include <mpi.h>

double **seidel(double border_x, double border_y,
                int Nx, int Ny,
                int local_Nx,
                double eps,
                double (*f)(double, double),
                double (*u_left)(double), double (*u_right)(double),
                double (*u_up)(double), double (*u_down)(double),
                int max_iter,
                int size, int rank);

#endif