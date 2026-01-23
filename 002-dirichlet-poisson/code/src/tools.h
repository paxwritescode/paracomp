#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>
#include <stdio.h>

double **alloc_matrix(int m, int n);
double max_in_matrix_diff(double** A, double** B, int m, int n);
void free_matrix(double **A, int n);
void init_matrix(double** A, int m, int n, double value);

#endif