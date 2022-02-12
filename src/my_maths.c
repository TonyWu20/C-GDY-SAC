#include "my_maths.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI (atan(1) * 4)

double *get_nth_vector_from_matrix(double *matrix, int rowSize, int colSize,
                                   int nth)
{
    double *vec = NULL;
    if (nth >= rowSize)
    {
        printf("Slice out of matrix range.\n");
        return vec;
    }
    vec = malloc(sizeof(double) * colSize);
    for (int i = 0; i < colSize; ++i)
    {
        vec[i] = matrix[nth * colSize + i];
    }
    return vec;
}

double *create_vector_a_to_b(double *a, double *b, int vectorSize)
{
    double *c = malloc(sizeof(double) * vectorSize);
    memcpy(c, b, sizeof(double) * vectorSize);
    cblas_daxpy(3, -1, a, 1, c, 1);
    return c;
}

void normalize_vector(double *u)
{
    double norm = cblas_dnrm2(3, u, 1);
    cblas_dscal(3, 1 / norm, u, 1);
}

double dot_product(double *u, double *v)
{
    double ret = cblas_ddot(3, u, 1, v, 1);
    return ret;
}

double vector_angle(double *u, double *v)
{
    double dotProduct = dot_product(u, v);
    double u_norm = cblas_dnrm2(3, u, 1);
    double v_norm = cblas_dnrm2(3, v, 1);
    double angle = acos(dotProduct / (u_norm * v_norm));
    return angle;
}

double *centroid_of_points(double *coords, int rowSize, int colSize)
{
    /* ColMajor */
    double avg_vec[colSize];
    for (int i = 0; i < colSize; ++i)
    {
        avg_vec[i] = (double)1 / (double)colSize;
    }
    double *centroid = calloc(4, sizeof(double));
    cblas_dgemv(CblasColMajor, CblasNoTrans, rowSize, colSize, 1, coords,
                rowSize, avg_vec, 1, 0, centroid, 1);
    return centroid;
}

double *cross_product(double *a, double *b) // Return normalized vector
{
    double a_skew_mat[16] = {0,     a[2],  -a[1], 0, //
                             -a[2], 0,     a[0],  0, //
                             a[1],  -a[0], 0,     0, 0, 0, 0, 1};
    double *cProduct = calloc(4, sizeof(double));
    cblas_dgemv(CblasColMajor, CblasNoTrans, 4, 4, 1, a_skew_mat, 4, b, 1, 0,
                cProduct, 1);
    return cProduct;
}

// Determine the transform matrix to make u parallel to v
double *rotate_u_to_v(double *u, double *v)
{
}

double *rotate_angle_around_axis(double *axis, double rad)
{
}

double *translate_mat_a_to_b(double *src, double *dest)
{
}

int roundupBiggerTenth(int number)
{
    int rounded = (number / 10 + 1) * 10;
    return rounded;
}
