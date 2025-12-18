#include "matrices.h"

double **alloc_matrix(int n)
{
    int i = 0;
    double** y = calloc(n, sizeof(*y));

    if (!y)
    {
        printf("Error: memory not allocated!\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < n; i++)
    {
        y[i] = calloc(n, sizeof(*y[i]));
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

void free_matrix(double** y, int n)
{
    int i = 0;
    for(i = 0; i < n; i++)
    {
        free(y[i]);
    }
    free(y);
}

double** substract_matrix(int n, double** y_1, double** y_2)
{
    int i = 0, j = 0;
    double** res = alloc_matrix(n);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            res[i][j] = y_1[i][j] - y_2[i][j];
        }
    }

    return res;
}

double compute_matrix_norm(int n, double** y)
{
    int i = 0, j = 0;
    double inf_norm = 0.0;


    for(i = 0; i < n; i++)
    {
        double cur_inf_norm = 0.0;
        for (j = 0; j < n; j++)
        {
            cur_inf_norm += fabs(y[i][j]);
        }
        if (cur_inf_norm > inf_norm)
        {
            inf_norm = cur_inf_norm;
        }
    }
    
    return inf_norm;
}