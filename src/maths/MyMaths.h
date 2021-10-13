#include "../main.h"
#include <math.h>
#include <stdlib.h>
void multiplyMatrices_mx3(double first[][3], double sec[][3],
                          double result[][3], int row1);
double VecAngle(double u[3], double v[3]);
double dotProduct(double u[], double v[], int size);
void initVector(double coord_a[3], double coord_b[3], double *result);
double NormVector(double u[], int size);
void moveMatrix(double Mat[][3], double *u, int ElmNum);
