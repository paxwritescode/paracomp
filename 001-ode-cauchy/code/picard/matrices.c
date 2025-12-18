#include "matrices.h"

double **alloc_matrix(int n, int m)
{
    int i = 0;
    double **y = calloc(n, sizeof(*y));

    if (!y)
    {
        printf("Error: memory not allocated!\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < n; i++)
    {
        y[i] = calloc(m, sizeof(*y[i]));
        if (!y[i])
        {
            for (int j = 0; j < i; j++)
            {
                free(y[j]);
            }
            free(y);
            exit(EXIT_FAILURE);
        }
    }
    return y;
}

void free_matrix(double **y, int n)
{
    int i = 0;
    for (i = 0; i < n; i++)
    {
        free(y[i]);
    }
    free(y);
}

double compute_diff_norm(int n, int m, double **y1, double **y2)
{
    double norm = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double row_sum = 0.0;
        for (int j = 0; j < m; ++j)
            row_sum += fabs(y1[i][j] - y2[i][j]);

        if (row_sum > norm)
            norm = row_sum;
    }
    return norm;
}