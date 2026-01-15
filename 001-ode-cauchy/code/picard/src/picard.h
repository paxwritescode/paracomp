#ifndef PICARD_H
#define PICARD_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double **picard_method(int n, double **f, double **A, double eps, double t_0, double t, int m, FILE *f_diff, FILE *f_norm);

#endif