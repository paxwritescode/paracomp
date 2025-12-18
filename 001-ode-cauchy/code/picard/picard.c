#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "picard.h"
#include "generate.h"
#include "matrices.h"

double **picard_method(int n, double *f, double eps, double t_0, double t, int m)
{
    double **y_prev = alloc_matrix(n, m);
    double **y = alloc_matrix(n, m);
    int i = 0, j = 0;
    int iterations = 0;
    double delta_t = (t - t_0) / (double)m;

    /* Iterations of Picard method */

    free_matrix(y_prev, n);
    return y;
}