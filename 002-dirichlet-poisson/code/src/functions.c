#include "functions.h"

double f(double x, double y)
{
    return (2 - y * y) * cos(x);
}

double u_0y(double y)
{
    return y * y;
}

double u_x0(double x)
{
    return 0;
}

double u_piy(double y)
{
    return - y * y;
}

double u_x1(double x)
{
    return cos(x);
}