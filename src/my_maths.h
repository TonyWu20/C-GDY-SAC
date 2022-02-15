#pragma once
#include <Accelerate/Accelerate.h>
#include <math.h>
#include <simd/simd.h>

#define PI (atan(1) * 4)

/* Get the ith vector from the mxn matrix as with the dimension of 1xn
 *
 * arguments:
 * double *matrix, int rowSize, int colSize, int nth (start from 0)
 *
 * returns:
 * A malloced double * with size of colSize
 */
double *get_nth_vector_from_matrix(double *matrix, int rowSize, int colSize,
                                   int nth);
simd_double4 get_simd_double4_from_matrix(double *matrix, int rowSize,
                                          int colSize, int i);

/* Create a vector point from a to b
 *
 * arguments:
 * double *a, double *b, int vectorSize
 *
 * returns:
 * A malloced double * with size of vectorSize
 */
double *create_vector_a_to_b(double *a, double *b, int vectorSize);

/* normalized vector in place */
void normalize_vector(double *vec);

/* Call cblas_ddot to calculate dot_product of the xyz part */
double dot_product(double *u, double *v);

/* Calculate angle of two vector */
double vector_angle(double *u, double *v);

double simd_vector_angle(simd_double3 u, simd_double3 v);

/* Calculate centroid of the given set of vertices */
double *centroid_of_points(double *coords, int rowSize, int colSize);

simd_double3 simd_centroid_of_points(simd_double3 *vecArray, int vecArraySize);
/* Return cross_product by matrix-vector multiplication
 *
 * returns:
 * A malloced double * 1x4 vector, not normalized
 */
double *cross_product(double *a, double *b);

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
