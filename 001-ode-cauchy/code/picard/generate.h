#include <math.h>

#ifndef GENERATE_H
#define GENERATE_H

void generate_test_rhs(double t, double *f);
void generate_rhs(int n, double t, double *f);
int generate_threads(int *threads, int max_threads);

#endif