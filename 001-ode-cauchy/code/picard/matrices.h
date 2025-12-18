#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MATRICES_H
#define MATRICES_H

double **alloc_matrix(int n);
void free_matrix(double** y, int n);
double** substract_matrix(int n, double** y_1, double** y_2);
double compute_matrix_norm(int n, double** y);

#endif