#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "picard.h"
#include "generate.h"
#include "matrices.h"


double **picard_method(int n, double *f, double eps, double t_0, double t)
{
    double** y_prev = alloc_matrix(n);
    double** y = alloc_matrix(n);

    /* Iterations of Picard method */
    
    free_matrix(y_prev, n);
    return y;
}