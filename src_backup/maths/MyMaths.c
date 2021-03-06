#include "MyMaths.h"
void multiplyMatrices_mx3(double first[][3], double sec[][3],
                          double result[][3], int row1)
{
    /* Initializing elements of matrix mult to 0 */
    for (int i = 0; i < row1; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = 0;
        }
    }

    /* Multiply first and second matrices and storing it into result */
    for (int i = 0; i < row1; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int n = 0; n < 3; n++)
            {
                result[i][j] += first[i][n] * sec[n][j];
            }
        }
    }
}
void initVector(double coord_a[3], double coord_b[3], double result[3])
{
    for (int i = 0; i < 3; i++)
    {
        result[i] = coord_b[i] - coord_a[i];
    }
}

void freeVector(double *u)
{
    if (u != NULL)
        free(u);
}

double dotProduct(double u[], double v[], int size)
{
    double result = 0.0;
    for (int i = 0; i < size; i++)
    {
        result += u[i] * v[i];
    }
    return result;
}

/** Calculate the angle between two 3d vector, return in degree;
 * Args:
 * 		double u[3];
 * 		double v[3];
 * returns:
 * 		double angle (deg);
 */
double VecAngle(double u[3], double v[3])
{
    double angle = 0.0;
    double dot_uv = 0.0;
    double norm_u, norm_v;
    norm_u = norm_v = 0.0;
    dot_uv = dotProduct(u, v, 3);
    norm_u = NormVector(u, 3);
    norm_v = NormVector(v, 3);
    angle = acos(dot_uv / (norm_u * norm_v));
    return angle * 180 / M_PI;
}

double NormVector(double u[], int size)
{
    double norm = 0.0;
    for (int i = 0; i < size; i++)
    {
        norm += u[i] * u[i];
    }
    norm = sqrt(norm);
    return norm;
}

/** Add Matrix by vector u
 *  Args:
 *  	double **Mat; arrays of coords,
 *  	double *u; vector,
 *  	int ElmNum; number of elements in Mat;
 */
void moveMatrix(double Mat[][3], double *u, int ElmNum)
{
    for (int i = 0; i < ElmNum; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Mat[i][j] += u[j];
        }
    }
}
