#include <math.h>

#ifndef GENERATE_H
#define GENERATE_H

double **generate_rhs(int n, int m, double t0, double t);
int generate_threads(int *threads, int max_threads);
double **generate_matrix(int n);

#endif