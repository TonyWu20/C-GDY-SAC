#include "my_maths.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI (atan(1) * 4)

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

// Determine the transform matrix to make u parallel to v
int roundupBiggerTenth(int number)
{
    int rounded = (number / 10 + 1) * 10;
    return rounded;
}
