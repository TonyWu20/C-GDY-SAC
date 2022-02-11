#include "my_maths.h"
#include <stdio.h>

double molCoords[] = {
    -0.04576804, 0.04703987,  0, 1, 0.95426804,  0.04456013, 0, 1,
    -0.61844704, -0.93819621, 0, 1, -0.61305652, 1.03510268, 0, 1,
    1.52155652,  -0.94350268, 0, 1, 1.52694704,  1.02979621, 0, 1,
};

/* Test retrieve vector, dot_product and vector_angle */
void test_1()
{
    double *a = get_nth_vector_from_matrix(molCoords, 6, 4, 0);
    double *b = get_nth_vector_from_matrix(molCoords, 6, 4, 1);
    double *c = get_nth_vector_from_matrix(molCoords, 6, 4, 2);
    double *ba = create_vector_a_to_b(a, b, 4);
    free(b);
    double *ca = create_vector_a_to_b(a, c, 4);
    free(a);
    free(c);
    double angle_ba_ca = vector_angle(ba, ca);
    printf("%f\n", angle_ba_ca * 180 / PI);
    free(ba);
    free(ca);
}

void test_centroid()
{
    double *centroid = centroid_of_points(molCoords, 4, 6);
    printf("%f %f %f %f\n", centroid[0], centroid[1], centroid[2], centroid[3]);
    free(centroid);
}

int main(int argc, char *argv[])
{
    test_centroid();
    return 0;
}
