#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double molCoords[] = {
    -0.04576804, 0.04703987,  0, 1, 0.95426804,  0.04456013, 0, 1,
    -0.61844704, -0.93819621, 0, 1, -0.61305652, 1.03510268, 0, 1,
    1.52155652,  -0.94350268, 0, 1, 1.52694704,  1.02979621, 0, 1,
};

/* Test retrieve vector, dot_product and vector_angle */
void test_simd()
{
    simd_double4 a = get_simd_double4_from_matrix(molCoords, 6, 4, 0);
    simd_double4 b = get_simd_double4_from_matrix(molCoords, 6, 4, 1);
    simd_double4 c = get_simd_double4_from_matrix(molCoords, 6, 4, 2);
    simd_double3 ba = simd_make_double3(b - a);
    simd_double3 ca = simd_make_double3(c - a);
    simd_double3 yAxis = {0, 1, 0};
    simd_double3 pNormal = simd_normalize(simd_cross(ba, ca));
    double angle = acos(simd_dot(pNormal, yAxis) / simd_length(pNormal) *
                        simd_length(yAxis));
    printf("%f\n", angle * 180 / PI);
    simd_quatd rotationQuat = simd_quaternion(angle, ba);
    simd_double3 vecArray[6];
    for (int i = 0; i < 6; ++i)
    {
        vecArray[i] =
            simd_make_double3(get_simd_double4_from_matrix(molCoords, 6, 4, i));
        vecArray[i] = simd_act(rotationQuat, vecArray[i]);
        printf("%f %f %f\n", vecArray[i].x, vecArray[i].y, vecArray[i].z);
    }
}
/* void test_1() */
/* { */
/*     double *a = get_nth_vector_from_matrix(molCoords, 6, 4, 0); */
/*     double *b = get_nth_vector_from_matrix(molCoords, 6, 4, 1); */
/*     double *c = get_nth_vector_from_matrix(molCoords, 6, 4, 2); */
/*     double *ba = create_vector_a_to_b(a, b, 4); */
/*     free(b); */
/*     double *ca = create_vector_a_to_b(a, c, 4); */
/*     free(a); */
/*     free(c); */
/*     double *cProduct_ba_ca = cross_product(ba, ca); */
/*     double yAxis[4] = {0, 1, 0, 1}; */
/*     double angle_ba_ca = vector_angle(cProduct_ba_ca, yAxis); */
/*     printf("%f\n", angle_ba_ca * 180 / PI); */
/*     double *rotated; */
/*     rotate_coords_around_axis_by_angle(ba, angle_ba_ca, molCoords, 4, 6, */
/*                                        &rotated); */
/*     for (int i = 0; i < 6; ++i) */
/*     { */
/*         for (int j = 0; j < 4; ++j) */
/*         { */
/*             printf("%f ", rotated[i * 4 + j]); */
/*         } */
/*         printf("\n"); */
/*     } */
/*     free(ba); */
/*     free(ca); */
/*     free(cProduct_ba_ca); */
/*     free(rotated); */
/* } */
/*  */
void test_centroid()
{
    clock_t t;
    t = clock();
    int i = 0;
    while (i < 100000)
    {
        double *centroid = centroid_of_points(molCoords, 4, 6);
        free(centroid);
        ++i;
    }
    t = clock() - t;
    double time_taken = ((double)t) / CLOCKS_PER_SEC;
    printf("BLAS centroid took %f seconds to execute\n", time_taken);
}

void test_simd_centroid()
{
    clock_t t;
    t = clock();
    int i = 0;
    simd_double3 vecArray[6];
    for (int i = 0; i < 6; ++i)
    {
        vecArray[i] =
            simd_make_double3(get_simd_double4_from_matrix(molCoords, 6, 4, i));
    }
    simd_double3 centroid;
    while (i < 100000)
    {
        centroid = simd_centroid_of_points(vecArray, 6);
        ++i;
    }
    t = clock() - t;
    double time_taken = ((double)t) / CLOCKS_PER_SEC;
    printf("%f %f %f\n", centroid[0], centroid[1], centroid[2]);
    printf("SIMD centroid took %f seconds to execute\n", time_taken);
}
int main(int argc, char *argv[])
{
    clock_t t;
    t = clock();
    test_centroid();
    t = clock() - t;
    double time_taken = ((double)t) / CLOCKS_PER_SEC;
    printf("test_centroid took %f seconds to execute\n", time_taken);
    t = clock();
    test_simd();
    t = clock() - t;
    time_taken = ((double)t) / CLOCKS_PER_SEC;
    printf("test_simd took %f seconds to execute\n", time_taken);
    test_simd_centroid();
    return 0;
}
