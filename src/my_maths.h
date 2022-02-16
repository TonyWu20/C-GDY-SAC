#pragma once
#include <math.h>
#include <simd/simd.h>

#define PI (atan(1) * 4)

simd_double4 get_simd_double4_from_matrix(double *matrix, int rowSize,
                                          int colSize, int i);

double simd_vector_angle(simd_double3 u, simd_double3 v);

/* Calculate centroid of the given set of vertices */

simd_double3 simd_centroid_of_points(simd_double3 *vecArray, int vecArraySize);

/* Matrix *translate_mat_a_to_b(double *center_a, double *center_b); */
/* void translate_a_to_b(Matrix *trans_mat, Matrix *coords, Matrix **result); */
/*  */
/* Matrix *fractionalCoordMatrix(Matrix *lat_vectors); */

simd_double3x3 fracCoordMat(simd_double3x3 latVectors);

simd_double4x4 translateMatrix(simd_double3 fromVec, simd_double3 toVec);

/* Round up number to the bigger nearest tenth.
 * E.g.: 374 -> 380; 376 -> 380;
 * Args: number, double type
 * Returns: rounded double
 */
int roundupBiggerTenth(int number);
