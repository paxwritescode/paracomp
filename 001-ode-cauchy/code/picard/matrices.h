#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MATRICES_H
#define MATRICES_H

double **alloc_matrix(int n, int m);
void free_matrix(double **y, int n);
double compute_diff_norm(int n, int m, double **y1, double **y2);
double* matrix_mul_vector(int n, double** A, double* y_m);

#endif