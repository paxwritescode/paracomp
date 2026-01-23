#include "tools.h"

double **alloc_matrix(int m, int n)
{
    double **V = calloc(m, sizeof(double*));
    if (!V)
    {
        printf("Error: memory not allocated!\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < m; i++)
    {
        V[i] = calloc(n, sizeof(double));
        if (!V[i])
        {
            printf("Error: memory not allocated!\n");
            exit(EXIT_FAILURE);
        }
    }
    return V;
}

double max_in_matrix_diff(double** A, double** B, int m, int n)
{
    double max = A[0][0] - B[0][0];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            if (fabs(A[i][j] - B[i][j]) > max)
                max = A[i][j];

    return max;
}

void free_matrix(double **A, int n)
{
    int i = 0;
    for (i = 0; i < n; i++)
    {
        free(A[i]);
    }
    free(A);
}

void init_matrix(double** A, int m, int n, double value)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = value;
}