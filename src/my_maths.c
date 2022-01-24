#include "my_maths.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI (atan(1) * 4)

Matrix *matrix_view_array(double base[], int m, int n)
{
    Matrix *ret = create_matrix(m, n);
    for (int i = 0; i < m * n; ++i)
    {
        ret->value[i / n][i % n] = base[i];
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
        res += m->value[i][0] * m->value[i][0];
    }
    norm = sqrt(res);
    return norm;
}

void normalize_vector(Matrix *vec)
{
    double norm = norm_of_vector(vec);
    multiply_matrix_with_scalar(vec, 1 / norm);
}

double dot_product(Matrix *u, Matrix *v)
{
    if (u->lines != v->lines)
    {
        printf("Inconsistent row sizes of u and v\n");
        return -1.0;
    }
    if (u->columns != v->columns)
    {
        printf("Inconsistent column sizes of u and v\n");
        return -1.0;
    }
    Matrix *v3_u, *v3_v;
    v3_u = create_matrix(3, 1);
    v3_v = create_matrix(3, 1);
    for (int i = 0; i < 3; ++i)
    {
        v3_u->value[i][0] = u->value[i][0];
        v3_v->value[i][0] = v->value[i][0];
    }
    Matrix *res;
    Matrix *u_t;
    get_transpose(v3_u, &u_t);
    multiply_matrices(u_t, v3_v, &res);
    double ret = res->value[0][0];
    destroy_matrix(res);
    destroy_matrix(u_t);
    destroy_matrix(v3_u);
    destroy_matrix(v3_v);
    free(res);
    free(u_t);
    free(v3_u);
    free(v3_v);
    return ret;
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
    double rot_x[] = {1, 0,       0,       0, 0, cos_rad, -sin_rad, 0,
                      0, sin_rad, cos_rad, 0, 0, 0,       0,        1};
    double rot_y[] = {cos_rad,  0, sin_rad, 0, 0, 1, 0, 0,
                      -sin_rad, 0, cos_rad, 0, 0, 0, 0, 1};
    double rot_z[] = {cos_rad, -sin_rad, 0, 0, sin_rad, cos_rad, 0, 0,
                      0,       0,        1, 0, 0,       0,       0, 1};
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
    double *ans = calloc(3, sizeof(double));
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

void rotate_around_origin(Matrix *rot_mat, Matrix *coords, Matrix **result)
{
    double *centroid = centroid_of_points(coords);
    for (int i = 0; i < 3; ++i)
    {
        rot_mat->value[i][3] = -centroid[i];
    }
    Matrix *tmp;
    multiply_matrices(rot_mat, coords, &tmp);
    double trans_m[] = {1, 0, 0, centroid[0], 0, 1, 0, centroid[1],
                        0, 0, 1, centroid[2], 0, 0, 0, 1};
    Matrix *trans_mat = matrix_view_array(trans_m, 4, 4);
    multiply_matrices(trans_mat, tmp, result);
    // Operations done. Tidy up memory
    free(centroid);
    destroy_matrix(tmp);
    destroy_matrix(trans_mat);
    free(tmp);
    free(trans_mat);
}

Matrix *cross_product(Matrix *a, Matrix *b) // Return normalized vector
{
    double a1 = a->value[0][0], a2 = a->value[1][0], a3 = a->value[2][0];
    double tmp_ca[] = {0,   -a3, a2,  0, //
                       a3,  0,   -a1, 0, //
                       -a2, a1,  0,   0, //
                       0,   0,   0,   1};
    Matrix *c_a = matrix_view_array(tmp_ca, 4, 4);
    Matrix *cross_product = NULL;
    multiply_matrices(c_a, b, &cross_product);
    normalize_vector(cross_product);
    cross_product->value[3][0] = 1;
    destroy_matrix(c_a);
    free(c_a);
    return cross_product;
}

// Determine the transform matrix to make u parallel to v
Matrix *rotate_u_to_v(Matrix *u, Matrix *v)
{
    Matrix *n = cross_product(u, v);
    double n1, n2, n3;
    n1 = n->value[0][0];
    n2 = n->value[1][0];
    n3 = n->value[2][0];
    double angle = vector_angle(u, v);
    double a = cos(angle);
    double b = 1 - a;
    double c = sin(angle);
    double Rij_n_T[] = {a + n1 * n1 * b,
                        n1 * n2 * b - n3 * c,
                        n1 * n3 * b + n2 * c,
                        0,
                        n1 * n2 * b + n3 * c,
                        a + n2 * n2 * b,
                        n2 * n3 * b - n1 * c,
                        0,
                        n1 * n3 * b - n2 * c,
                        n2 * n3 * b + n1 * c,
                        a + n3 * n3 * b,
                        0, //
                        0,
                        0,
                        0,
                        1};
    Matrix *T = matrix_view_array(Rij_n_T, 4, 4);
    destroy_matrix(n);
    free(n);
    return T;
}

Matrix *translate_mat_a_to_b(double *src, double *dest)
{
    Matrix *translate;
    double tx, ty, tz;
    tx = dest[0] - src[0];
    ty = dest[1] - src[1];
    tz = dest[2] - src[2] + 1.5;
    double t_array[] = {1, 0, 0, tx, 0, 1, 0, ty, 0, 0, 1, tz, 0, 0, 0, 1};
    translate = matrix_view_array(t_array, 4, 4);
    return translate;
}

void translate_a_to_b(Matrix *trans_mat, Matrix *coords, Matrix **result)
{
    multiply_matrices(trans_mat, coords, result);
}
