#include "my_maths.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

#define PI (atan(1) * 4)

double norm_of_vector(Matrix *m)
{
    double norm;
    double res = 0;
    for (int i = 0; i < 3; ++i)
    {
        res += pow(m->value[i][0], 2);
    }
    norm = sqrt(res);
    return norm;
}

double dot_product(Matrix *u, Matrix *v)
{
    if (u->columns != v->columns)
    {
        printf("Inconsistent column size of u and v");
        return -1.0;
    }
    double res = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        res += u->value[i][0] * v->value[i][0];
    }
    return res;
}

double vector_angle(Matrix *u, Matrix *v)
{
    double dot_uv = dot_product(u, v);
    double norm_u = norm_of_vector(u);
    double norm_v = norm_of_vector(v);
    double angle = acos(dot_uv / (norm_u * norm_v));
    return angle;
}

Matrix *rotationMatrix(double rad, char axis)
{
    double cos_rad = cos(rad);
    double sin_rad = sin(rad);
    Matrix *rotMat = create_matrix(4, 4);
    double rot_x[][4] = {{1, 0, 0, 0},
                         {0, cos_rad, -1 * sin_rad, 0},
                         {0, sin_rad, cos_rad, 0},
                         {0, 0, 0, 1}};
    double rot_y[][4] = {{cos_rad, 0, sin_rad, 0},
                         {0, 1, 0, 0},
                         {-1 * sin_rad, 0, cos_rad, 0},
                         {0, 0, 0, 1}};
    double rot_z[][4] = {{cos_rad, -1 * sin_rad, 0, 0},
                         {sin_rad, cos_rad, 0, 0},
                         {0, 0, 1, 0},
                         {0, 0, 0, 1}};
    switch (axis)
    {
    case 'X':
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                rotMat->value[i][j] = rot_x[i][j];
            }
        }
        break;
    }
    case 'Y':
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                rotMat->value[i][j] = rot_y[i][j];
            }
        }
        break;
    }
    case 'Z':
    {
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                rotMat->value[i][j] = rot_z[i][j];
            }
        }
        break;
    }
    default:
        break;
    }
    return rotMat;
}

double cross_product(Matrix *u, Matrix *v)
{
}
