#pragma once
#include "matrix.h"
#include <math.h>

#define PI (atan(1) * 4)
Matrix *matrix_view_array(double base[][4], int m, int n);
Matrix *col_vector_view_array(double base[], int m);
double norm_of_vector(Matrix *m);
double dot_product(Matrix *u, Matrix *v);
double vector_angle(Matrix *u, Matrix *v);
Matrix *cross_product(Matrix *a, Matrix *b);
Matrix *rotationMatrix(double rad, char axis);
double *centroid_of_points(Matrix *coords);
void rotate_around_origin(Matrix *coords, double rad, char axis,
                          Matrix **result);
Matrix *cross_product(Matrix *a, Matrix *b); // Return normalized vector
