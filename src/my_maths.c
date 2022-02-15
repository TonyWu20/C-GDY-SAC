#include "my_maths.h"
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

double dot_product(double *u, double *v)
{
    return cblas_ddot(3, u, 1, v, 1);
}

double vector_angle(double *u, double *v)
{
    double dP = dot_product(u, v);
    double u_norm = sqrt(dot_product(u, u));
    double v_norm = sqrt(dot_product(v, v));
    return acos(dP / (u_norm * v_norm));
}

double simd_vector_angle(simd_double3 u, simd_double3 v)
{
    double dP_uv = simd_dot(u, v);
    return acos(dP_uv / (simd_length(u) * simd_length(v)));
}

simd_double4 get_simd_double4_from_matrix(double *matrix, int rowSize,
                                          int colSize, int i)
{
    simd_double4 vec = simd_double((simd_double4){0, 0, 0, 0});
    if (i >= rowSize)
    {
        printf("Slice out of matrix range.\n");
        return vec;
    }
    vec = simd_double(
        (simd_double4){matrix[i * colSize + 0], matrix[i * colSize + 1],
                       matrix[i * colSize + 2], matrix[i * colSize + 3]});
    return vec;
}

simd_double3 simd_centroid_of_points(simd_double3 *vecArray, int vecArraySize)
{
    simd_double3 res = {0, 0, 0};
    for (int i = 0; i < vecArraySize; ++i)
    {
        res += vecArray[i];
    }
    res = res / (double)vecArraySize;
    return res;
}

simd_double3x3 fracCoordMat(simd_double3x3 latVectors)
{
    simd_double3 a = latVectors.columns[0];
    simd_double3 b = latVectors.columns[1];
    simd_double3 c = latVectors.columns[2];
    double aLen = simd_length(a);
    double bLen = simd_length(b);
    double cLen = simd_length(c);
    double alpha = simd_vector_angle(b, c);
    double beta = simd_vector_angle(a, c);
    double gamma = simd_vector_angle(a, b);
    double vol = simd_dot(a, simd_cross(b, c));
    simd_double3 x_cart = {aLen, 0, 0};
    simd_double3 y_cart = {bLen * cos(gamma), bLen * sin(gamma), 0};
    simd_double3 z_cart = {cLen * cos(beta),
                           cLen * (cos(alpha) - cos(beta) * cos(gamma)) /
                               sin(gamma),
                           vol / (aLen * bLen * sin(gamma))};
    simd_double3x3 toCart = simd_matrix(x_cart, y_cart, z_cart);
    simd_double3x3 toFrac = simd_inverse(toCart);
    return toFrac;
}

simd_double4x4 translateMatrix(simd_double3 fromVec, simd_double3 toVec)
{
    simd_double3 translation = toVec - fromVec;
    simd_double4 diagonal = {1, 1, 1, 1};
    simd_double4x4 transMat = simd_diagonal_matrix(diagonal);
    transMat.columns[3].x = translation.x;
    transMat.columns[3].y = translation.y;
    transMat.columns[3].z = translation.z;
    return transMat;
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

double *centroid_of_points(double *coords, int rowSize, int colSize)
{
    double avg_vec[colSize];
    for (int i = 0; i < colSize; ++i)
    {
        avg_vec[i] = 1 / (double)colSize;
    }
    double *centroid = calloc(4, sizeof(double));
    cblas_dgemv(CblasColMajor, CblasNoTrans, rowSize, colSize, 1, coords,
                rowSize, avg_vec, 1, 0, centroid, 1);
    return centroid;
}

void normalize_vector(double *vec)
{
    double norm = sqrt(cblas_ddot(3, vec, 1, vec, 1));
    cblas_dscal(3, 1 / norm, vec, 1);
}

// Determine the transform matrix to make u parallel to v

void rotate_coords_around_axis_by_angle(double *axis, double rad,
                                        double *coords, int rowSize,
                                        int colSize, double **result)
{
    normalize_vector(axis);
    double ux, uy, uz;
    ux = axis[0];
    uy = axis[1];
    uz = axis[2];
    double a = cos(rad);
    double b = 1 - a;
    double c = sin(rad);
    double R_axis_rad[] = {a + ux * ux * b,
                           uy * ux * b + uz * c,
                           uz * ux * b + uy * c,
                           0, //
                           ux * uy * b - uz * c,
                           a + uy * uy * b,
                           uz * uy * b + ux * c,
                           0, //
                           ux * uz * b + uy * c,
                           uy * uz * b - ux * c,
                           a + uz * uz * b,
                           0,
                           0,
                           0,
                           0,
                           1};
    *result = calloc(rowSize * colSize, sizeof(double));
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rowSize, colSize, 4,
                1, R_axis_rad, 4, coords, 4, 0, *result, 4);
    /* return result; */
}

int roundupBiggerTenth(int number)
{
    int rounded = (number / 10 + 1) * 10;
    return rounded;
}
