#include "generate.h"

void generate_test_rhs(double t, double *f)
{
    f[0] = t;
    f[1] = t * t;
    f[2] = sin(t);
}

void generate_rhs(int n, double t, double *f)
{
    for (int i = 0; i < n; i++)
    {
        f[i] = sin(t + (double)i);
    }
}

int generate_threads(int *threads, int max_threads)
{
    int k = 0;
    int p = 1;

    while (p <= max_threads)
    {
        threads[k++] = p;
        p *= 2;
    }
    return k;
}