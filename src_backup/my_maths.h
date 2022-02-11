#pragma once
#include "matrix.h"
#include <math.h>

#define PI (atan(1) * 4)
Matrix *matrix_view_array(double base[], int m, int n);
Matrix *col_vector_view_array(double base[], int m);
double norm_of_vector(Matrix *m);
void normalize_vector(Matrix *vec);
double dot_product(Matrix *u, Matrix *v);
double vector_angle(Matrix *u, Matrix *v);
Matrix *cross_product(Matrix *a, Matrix *b);
Matrix *rotationMatrix(double rad, char axis);
double *centroid_of_points(Matrix *coords);
Matrix *cross_product(Matrix *a, Matrix *b); // Return normalized vector
Matrix *rotate_u_to_v(Matrix *u, Matrix *v);
Matrix *rotate_angle_around_axis(Matrix *axis, double rad);
Matrix *translate_mat_a_to_b(double *center_a, double *center_b);
void translate_a_to_b(Matrix *trans_mat, Matrix *coords, Matrix **result);

Matrix *fractionalCoordMatrix(Matrix *lat_vectors);

/* Round up number to the bigger nearest tenth.
 * E.g.: 374 -> 380; 376 -> 380;
 * Args: number, double type
 * Returns: rounded double
 */
int roundupBiggerTenth(int number);
