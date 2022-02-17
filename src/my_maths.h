#pragma once
#include <math.h>
#define PI (atan(1) * 4)
#define VEC_OVERLOAD __attribute__((__overloadable__))
#define VEC_INLINE __attribute__((__always_inline__))
#define VEC_NODEBUG __attribute__((__nodebug__))
#define VEC_CFUNC VEC_OVERLOAD VEC_INLINE VEC_NODEBUG
/* typedef union */
/* { */
/*     double xyzw[4]; */
/*     struct */
/*     { */
/*         double x; */
/*         double y; */
/*         double z; */
/*         double w; */
/*     }; */
/* } Vector_double4; */

/* typedef union */
/* { */
/*     double xyz[3]; */
/*     struct */
/*     { */
/*         double x, y, z; */
/*     }; */
/* } Vector_double3; */

typedef __attribute__((__ext_vector_type__(3),
                       __aligned__(16))) double Vector_double3;

/* typedef union */
/* { */
/*     Vector_double3 columns[3]; */
/*     struct */
/*     { */
/*         Vector_double3 i, j, k; */
/*     }; */
/* } Matrix_3x3; */
/*  */
/* typedef union */
/* { */
/*     Vector_double4 columns[4]; */
/*     struct */
/*     { */
/*         Vector_double4 i, j, k, e; */
/*     }; */
/* } Matrix_4x4; */
/* static inline VEC_CFUNC Vector_double3 make_Vector_double3(double x, double
 * y, */
/*                                                            double z) */
/* { */
/*     Vector_double3 result; */
/*     result.x = x; */
/*     result.y = y; */
/*     result.z = z; */
/*     return result; */
/* }; */
/* static inline VEC_CFUNC Vector_double3 make_Vector_double3(double vec[3]) */
/* { */
/*     Vector_double3 result; */
/*     result.x = vec[0]; */
/*     result.y = vec[1]; */
/*     result.z = vec[2]; */
/*     return result; */
/* } */
/* Matrix *matrix_view_array(double base[], int m, int n); */
/* Matrix *col_vector_view_array(double base[], int m); */
/* double norm_of_vector(Matrix *m); */
/* void normalize_vector(Matrix *vec); */
/* double dot_product(Matrix *u, Matrix *v); */
/* double vector_angle(Matrix *u, Matrix *v); */
/* Matrix *cross_product(Matrix *a, Matrix *b); */
/* Matrix *rotationMatrix(double rad, char axis); */
/* double *centroid_of_points(Matrix *coords); */
/* Matrix *cross_product(Matrix *a, Matrix *b); // Return normalized vector */
/* Matrix *rotate_u_to_v(Matrix *u, Matrix *v); */
/* Matrix *rotate_angle_around_axis(Matrix *axis, double rad); */
/* Matrix *translate_mat_a_to_b(double *center_a, double *center_b); */
/* void translate_a_to_b(Matrix *trans_mat, Matrix *coords, Matrix **result); */
/*  */
/* Matrix *fractionalCoordMatrix(Matrix *lat_vectors); */

/* Round up number to the bigger nearest tenth.
 * E.g.: 374 -> 380; 376 -> 380;
 * Args: number, double type
 * Returns: rounded double
 */
int roundupBiggerTenth(int number);
