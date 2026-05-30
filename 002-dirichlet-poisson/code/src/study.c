#include "study.h"

void study_time_vs_grid_size(int rank, int size)
{
    int test_grids[][2] = {
        {40, 10}, {40, 20}, {40, 40}, {40, 80}, {40, 100}, {40, 120}, {40, 140}, {40, 160}, {40, 180}, {40, 200}, 
        {100, 10}, {100, 20}, {100, 40}, {100, 80}, {100, 100}, {100, 120}, {100, 140}, {100, 160}, {100, 180}, {100, 200}, 
        {10, 40}, {20, 40}, {60, 40}, {80, 40}, {100, 40}, {120, 40}, {140, 40}, {160, 40}, {180, 40}, {200, 40}, 
        {10, 100}, {20, 100}, {60, 100}, {80, 100}, {100, 100}, {120, 100}, {140, 100}, {160, 100}, {180, 100}, {200, 100}, 
    };
    int num_tests = sizeof(test_grids) / sizeof(test_grids[0]);
    
    double eps = 1e-6;
    int max_iter = 200000;
    double right_border = PI, upper_border = 1.0;

    if (rank == 0 && size == 1)
    {
        FILE *f_out = fopen("results/time_vs_grid.csv", "w");
        if (f_out)
        {
            fprintf(f_out, "Nx,Ny,NumProcs,TimeSeconds\n");
            fclose(f_out);
        }
    }

    for (int t = 0; t < num_tests; t++)
    {
        int Nx = test_grids[t][0];
        int Ny = test_grids[t][1];

        MPI_Barrier(MPI_COMM_WORLD);
        double t_start = MPI_Wtime();

        double **V = seidel(right_border, upper_border, 
                            Nx, Ny, eps, 
                            &f, &u_0y, &u_piy, &u_x1, &u_x0, 
                            max_iter);

        double t_end = MPI_Wtime();
        double local_time = t_end - t_start;
        double global_time;

        MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            if (V) free_matrix(V, Nx + 1);
            
            FILE *f_out = fopen("results/time_vs_grid.csv", "a");
            if (f_out)
            {
                fprintf(f_out, "%d,%d,%d,%lf\n", Nx, Ny, size, global_time);
                fclose(f_out);
            }
        }
    }
}

void study_scaling_performance(int rank, int size)
{
    int Nx = 80;
    int Ny = 400;
    double eps = 1e-6;
    int max_iter = 100000;
    double right_border = PI, upper_border = 1.0;

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    double **V = seidel(right_border, upper_border, 
                        Nx, Ny, eps, 
                        &f, &u_0y, &u_piy, &u_x1, &u_x0, 
                        max_iter);

    double t_end = MPI_Wtime();
    double local_time = t_end - t_start;
    double global_time;

    MPI_Reduce(&local_time, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        if (V) free_matrix(V, Nx + 1);

        double t1 = -1.0;

        FILE *f_check = fopen("results/scaling_data.csv", "r");
        if (f_check)
        {
            char header[200];
            if (fgets(header, sizeof(header), f_check))
            {
                int p; double t, s, e;
                while (fscanf(f_check, "%d,%lf,%lf,%lf", &p, &t, &s, &e) == 4)
                {
                    if (p == 1)
                    {
                        t1 = t;
                        break;
                    }
                }
            }
            fclose(f_check);
        }

        if (size == 1)
        {
            t1 = global_time;
            FILE *f_out = fopen("results/scaling_data.csv", "w");
            if (f_out)
            {
                fprintf(f_out, "NumProcs,TimeSeconds,Speedup,Efficiency\n");
                fprintf(f_out, "%d,%lf,%lf,%lf\n", size, global_time, 1.0, 1.0);
                fclose(f_out);
            }
        }
        else
        {
            FILE *f_out = fopen("results/scaling_data.csv", "a");
            if (f_out)
            {
                double speedup = (t1 > 0.0) ? (t1 / global_time) : 0.0;
                double efficiency = (t1 > 0.0) ? (speedup / (double)size) : 0.0;
                
                fprintf(f_out, "%d,%lf,%lf,%lf\n", size, global_time, speedup, efficiency);
                fclose(f_out);
            }
        }
    }
}