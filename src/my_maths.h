#pragma once
#include <math.h>
#define PI (atan(1) * 4)
#define VEC_OVERLOAD __attribute__((__overloadable__))
#define VEC_INLINE __attribute__((__always_inline__))
#define VEC_NODEBUG __attribute__((__nodebug__))
#define VEC_CFUNC VEC_OVERLOAD VEC_INLINE VEC_NODEBUG

/* Definition */
typedef union
{
    double xyzw[4];
    struct
    {
        double x;
        double y;
        double z;
        double w;
    };
} vec_double4;

typedef union
{
    double xyz[3];
    struct
    {
        double x, y, z;
    };
} vec_double3;

typedef union
{
    vec_double3 columns[3];
    struct
    {
        vec_double3 i, j, k;
    };
} matrix_double3x3;

typedef union
{
    vec_double4 columns[4];
    struct
    {
        vec_double4 i, j, k, t;
    };
} matrix_double4x4;

typedef union
{
    double xyzr[4];
    struct
    {
        double x, y, z, r;
    };
} vec_quatd;

static inline VEC_CFUNC vec_double3 vec_make_double3(double x, double y,
                                                     double z);
static inline VEC_CFUNC vec_double3 vec_make_double3(double vec[3]);
static inline VEC_CFUNC vec_double3 vec_make_double3(vec_double4 other);
static inline VEC_CFUNC vec_double4 vec_make_double4(double x, double y,
                                                     double z, double w);
static inline VEC_CFUNC vec_double4 vec_make_double4(vec_double3 other,
                                                     double w);
static inline VEC_CFUNC vec_quatd vec_make_quaternion(double angle,
                                                      vec_double3 axis);
static inline VEC_CFUNC vec_quatd vec_make_quaternion(vec_double3 imag,
                                                      double real);
static inline VEC_CFUNC vec_quatd vec_make_quaternion(vec_double4 other);

static inline VEC_CFUNC double vec_dot(vec_double3 u, vec_double3 v);
static inline VEC_CFUNC double vec_length(vec_double3 u);
static inline VEC_CFUNC vec_double3 vec_normalize(vec_double3 u);
static inline VEC_CFUNC double vec_angle_uv(vec_double3 u, vec_double3 v);
static inline VEC_CFUNC vec_double3 vec_cross(vec_double3 u, vec_double3 v);
static inline VEC_CFUNC vec_double3 vec_act(vec_quatd q, vec_double3 v);
static inline VEC_CFUNC vec_double3 vec_imag(vec_quatd q);
static inline VEC_CFUNC double vec_real(vec_quatd q);
static inline VEC_CFUNC matrix_double4x4 vec_translate(vec_double3 from,
                                                       vec_double3 to);

/* Scale the vector
 * c = a*v;
 */
static inline VEC_CFUNC vec_double3 vec_mul(double a, vec_double3 v);
static inline VEC_CFUNC vec_double4 vec_mul(double a, vec_double4 v);
/* c = a*x + b*y */
static inline VEC_CFUNC vec_double3 vec_linear_comb(double a, vec_double3 x,
                                                    double b, vec_double3 y);
static inline VEC_CFUNC vec_double4 vec_linear_comb(double a, vec_double4 x,
                                                    double b, vec_double4 y);
/* c = x + y by vec_linear_comb */
static inline VEC_CFUNC vec_double3 vec_add(vec_double3 x, vec_double3 y);
static inline VEC_CFUNC vec_double4 vec_add(vec_double4 x, vec_double4 y);
/* c = x - y by vec_linear_comb */
static inline VEC_CFUNC vec_double3 vec_sub(vec_double3 x, vec_double3 y);
static inline VEC_CFUNC vec_double4 vec_sub(vec_double4 x, vec_double4 y);
/* Matrix vector multiplication */
static inline VEC_CFUNC vec_double3 vec_mul(matrix_double3x3 mat_3x3,
                                            vec_double3 u3);
static inline VEC_CFUNC vec_double4 vec_mul(matrix_double4x4 mat_4x4,
                                            vec_double4 u4);

/* Implementation */
static inline VEC_CFUNC vec_double3 vec_make_double3(double x, double y,
                                                     double z)
{
    vec_double3 result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
}
static inline VEC_CFUNC vec_double3 vec_make_double3(double vec[3])
{
    vec_double3 result;
    result.x = vec[0];
    result.y = vec[1];
    result.z = vec[2];
    return result;
}
static inline VEC_CFUNC vec_double3 vec_make_double3(vec_double4 other)
{
    vec_double3 result;
    for (int i = 0; i < 3; ++i)
    {
        result.xyz[i] = other.xyzw[i];
    }
    return result;
}
static inline VEC_CFUNC vec_double4 vec_make_double4(double x, double y,
                                                     double z, double w)
{
    vec_double4 result;
    result.x = x;
    result.y = y;
    result.z = z;
    result.w = w;
    return result;
}
static inline VEC_CFUNC vec_double4 vec_make_double4(vec_double3 other,
                                                     double w)
{
    vec_double4 result;
    for (int i = 0; i < 3; ++i)
    {
        result.xyzw[i] = other.xyz[i];
    }
    result.w = w;
    return result;
}
static inline VEC_CFUNC vec_quatd vec_make_quaternion(vec_double3 imag,
                                                      double real)
{
    vec_quatd result;
    for (int i = 0; i < 3; ++i)
    {
        result.xyzr[i] = imag.xyz[i];
    }
    result.r = real;
    return result;
}
static inline VEC_CFUNC vec_quatd vec_make_quaternion(vec_double4 other)
{
    vec_quatd result;
    for (int i = 0; i < 4; ++i)
    {
        result.xyzr[i] = other.xyzw[i];
    }
    return result;
}
static inline VEC_CFUNC vec_quatd vec_make_quaternion(double angle,
                                                      vec_double3 axis)
{
    vec_double3 imag;
    for (int i = 0; i < 3; ++i)
    {
        imag.xyz[i] = sin(angle / 2) * axis.xyz[i];
    }
    double real = cos(angle / 2);
    vec_quatd result = vec_make_quaternion(imag, real);
    return result;
}
static inline VEC_CFUNC vec_double3 vec_imag(vec_quatd q)
{
    vec_double3 imag;
    for (int i = 0; i < 3; ++i)
    {
        imag.xyz[i] = q.xyzr[i];
    }
    return imag;
}
static inline VEC_CFUNC double vec_real(vec_quatd q)
{
    return q.r;
}
static inline VEC_CFUNC vec_double3 vec_mul(double a, vec_double3 x)
{
    vec_double3 result;
    for (int i = 0; i < 3; ++i)
    {
        result.xyz[i] = x.xyz[i] * a;
    }
    return result;
}
static inline VEC_CFUNC vec_double4 vec_mul(double a, vec_double4 x)
{
    vec_double4 result;
    for (int i = 0; i < 4; ++i)
    {
        result.xyzw[i] = x.xyzw[i] * a;
    }
    return result;
}
static inline VEC_CFUNC vec_double3 vec_linear_comb(double a, vec_double3 x,
                                                    double b, vec_double3 y)
{
    vec_double3 result;
    vec_double3 aX = vec_mul(a, x), bY = vec_mul(b, y);
    for (int i = 0; i < 3; ++i)
    {
        result.xyz[i] = aX.xyz[i] + bY.xyz[i];
    }
    return result;
}
static inline VEC_CFUNC vec_double4 vec_linear_comb(double a, vec_double4 x,
                                                    double b, vec_double4 y)
{
    vec_double4 result;
    vec_double4 aX = vec_mul(a, x), bY = vec_mul(b, y);
    for (int i = 0; i < 4; ++i)
    {
        result.xyzw[i] = aX.xyzw[i] + bY.xyzw[i];
    }
    return result;
}
static inline VEC_CFUNC vec_double3 vec_add(vec_double3 x, vec_double3 y)
{
    return vec_linear_comb(1, x, 1, y);
}
static inline VEC_CFUNC vec_double4 vec_add(vec_double4 x, vec_double4 y)
{
    return vec_linear_comb(1, x, 1, y);
}
static inline VEC_CFUNC vec_double3 vec_sub(vec_double3 x, vec_double3 y)
{
    return vec_linear_comb(1, x, -1, y);
}
static inline VEC_CFUNC vec_double4 vec_sub(vec_double4 x, vec_double4 y)
{
    return vec_linear_comb(1, x, -1, y);
}
static inline VEC_CFUNC vec_double3 vec_mul(matrix_double3x3 mat_3x3,
                                            vec_double3 u3)
{
    vec_double3 result;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            result.xyz[i] = mat_3x3.columns[j].xyz[i] * u3.xyz[j];
        }
    }
    return result;
}
static inline VEC_CFUNC vec_double4 vec_mul(matrix_double4x4 mat_4x4,
                                            vec_double4 u4)
{
    vec_double4 result;
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 4; ++i)
        {
            result.xyzw[j] = mat_4x4.columns[i].xyzw[j] * u4.xyzw[i];
        }
    }
    return result;
}
static inline VEC_CFUNC matrix_double4x4 vec_translate(vec_double3 from,
                                                       vec_double3 to)
{
    vec_double3 fromTo = vec_sub(to, from);
    matrix_double4x4 transMat;
    transMat.i = vec_make_double4(1, 0, 0, 0);
    transMat.j = vec_make_double4(0, 1, 0, 0);
    transMat.k = vec_make_double4(0, 0, 1, 0);
    transMat.t = vec_make_double4(fromTo.x, fromTo.y, fromTo.z, 1);
    return transMat;
}
static inline VEC_CFUNC double vec_dot(vec_double3 u, vec_double3 v)
{
    double result = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        result += u.xyz[i] * v.xyz[i];
    }
    return result;
}

static inline VEC_CFUNC double vec_length(vec_double3 u)
{
    double result = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        result += u.xyz[i] * u.xyz[i];
    }
    return sqrt(result);
}

static inline VEC_CFUNC vec_double3 vec_normalize(vec_double3 u)
{
    double uLength = vec_length(u);
    vec_double3 result;
    for (int i = 0; i < 3; ++i)
    {
        result.xyz[i] = u.xyz[i] / uLength;
    }
    return result;
}
static inline VEC_CFUNC double vec_angle_uv(vec_double3 u, vec_double3 v)
{
    double dotProduct = vec_dot(u, v);
    double result = acos(dotProduct / (vec_length(u) * vec_length(v)));
    return result;
}
static inline VEC_CFUNC vec_double3 vec_cross(vec_double3 u, vec_double3 v)
{
    vec_double3 cross;
    cross.x = u.y * v.z - u.z * v.y; // a2b3 - a3b2
    cross.y = u.z * v.x - u.x * v.z; // a3b1 - a1b3
    cross.z = u.x * v.y - u.y * v.x; // a1b2 - a2b1
    return cross;
}
static inline VEC_CFUNC vec_double3 vec_centroid(vec_double3 *points,
                                                 int pointNums)
{
    vec_double3 result = {{0, 0, 0}};
    for (int i = 0; i < pointNums; ++i)
    {
        result.x += points[i].x;
        result.y += points[i].y;
        result.z += points[i].z;
    }
    result.x /= (double)pointNums;
    result.y /= (double)pointNums;
    result.z /= (double)pointNums;
    return result;
}
static inline VEC_CFUNC vec_double3 vec_act(vec_quatd q, vec_double3 v)
{
    vec_double3 cross = vec_cross(vec_imag(q), v);
    vec_double3 result;
    for (int i = 0; i < 3; ++i)
    {
        result.xyz[i] =
            v.xyz[i] + vec_real(q) * 2 * cross.xyz[i] + cross.xyz[i];
    }
    return result;
}
/* Matrix *translate_mat_a_to_b(double *center_a, double *center_b); */

/* Round up number to the bigger nearest tenth.

 * E.g.: 374 -> 380; 376 -> 380;

 * Args: number, double type

 * Returns: rounded double

 */
static inline int roundupBiggerTenth(int number);
static inline int roundupBiggerTenth(int number)
{
    int rounded = (number / 10 + 1) * 10;
    return rounded;
}
