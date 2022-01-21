#pragma once
#include "matrix.h"

double norm_of_vector(Matrix *m);
double dot_product(Matrix *u, Matrix *v);
double vector_angle(Matrix *u, Matrix *v);
Matrix *cross_product(Matrix *a, Matrix *b);
Matrix *rotationMatrix(double rad, char axis);
Matrix *center_of_mass(Matrix **vertices);
