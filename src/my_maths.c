#include "my_maths.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI (atan(1) * 4)

Matrix *matrix_view_array(double base[][4], int m, int n)
{
    Matrix *ret = create_matrix(m, n);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ret->value[i][j] = base[i][j];
        }
    }
    return ret;
}

Matrix *col_vector_view_array(double base[], int m)
{
    Matrix *ret = create_matrix(m, 1);
    for (int i = 0; i < m; ++i)
    {
        ret->value[i][0] = base[i];
    }
    return ret;
}

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
    Matrix *rotMat;
    switch (axis)
    {
    case 'X':
    {
        rotMat = matrix_view_array(rot_x, 4, 4);
        break;
    }
    case 'Y':
        rotMat = matrix_view_array(rot_y, 4, 4);
        break;
    case 'Z':
        rotMat = matrix_view_array(rot_z, 4, 4);
        break;
    default:
        rotMat = NULL;
        printf("Wrong axis value\n");
        break;
    }
    return rotMat;
}

double *centroid_of_points(Matrix *coords)
{
    double *ans = malloc(sizeof(double) * 3);
    for (int i = 0; i < coords->columns; ++i)
    {
        ans[0] += coords->value[0][i];
        ans[1] += coords->value[1][i];
        ans[2] += coords->value[2][i];
    }
    ans[0] /= coords->columns;
    ans[1] /= coords->columns;
    ans[2] /= coords->columns;
    return ans;
}

void rotate_around_origin(Matrix *coords, double rad, char axis,
                          Matrix **result)
{
    Matrix *rot_mat = rotationMatrix(rad, axis);
    double *centroid = centroid_of_points(coords);
    for (int i = 0; i < 3; ++i)
    {
        rot_mat->value[i][3] = -centroid[i];
    }
    Matrix *tmp;
    multiply_matrices(rot_mat, coords, &tmp);
    double trans_m[][4] = {{1, 0, 0, centroid[0]},
                           {0, 1, 0, centroid[1]},
                           {0, 0, 1, centroid[2]},
                           {0, 0, 0, 1}};
    Matrix *trans_mat = matrix_view_array(trans_m, 4, 4);
    multiply_matrices(trans_mat, tmp, result);
    // Operations done. Tidy up memory
    free(centroid);
    destroy_matrix(rot_mat);
    destroy_matrix(tmp);
    destroy_matrix(trans_mat);
    free(rot_mat);
    free(tmp);
    free(trans_mat);
}

Matrix *cross_product(Matrix *a, Matrix *b) // Return normalized vector
{
    double a1 = a->value[0][0], a2 = a->value[1][0], a3 = a->value[2][0];
    double tmp_ca[][4] = {
        {0, -a3, a2, 0}, {a3, 0, -a1, 0}, {-a2, a1, 0, 0}, {0, 0, 0, 1}};
    Matrix *c_a = matrix_view_array(tmp_ca, 4, 4);
    Matrix *cross_product = NULL;
    multiply_matrices(c_a, b, &cross_product);
    double norm_cross = norm_of_vector(cross_product);
    multiply_matrix_with_scalar(cross_product, 1 / norm_cross);
    destroy_matrix(c_a);
    free(c_a);
    return cross_product;
}
