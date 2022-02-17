#include "my_maths/my_maths.h"

// clang-format off
typedef double vec_double1;
typedef __attribute__((__ext_vector_type__(2))) double vec_double2;
typedef __attribute__((__ext_vector_type__(3), __aligned__(16))) double vec_double3;
typedef __attribute__((__ext_vector_type__(4), __aligned__(16))) double vec_double4;
/*! @abstract A matrix with 3 rows and 3 columns.*/
typedef struct { vec_double3 columns[3]; } vec_double3x3;
typedef struct { vec_double4 columns[4]; } vec_double4x4;
typedef struct { vec_double4 vector; } vec_quatd;
// clang-format on

static inline VEC_CFUNC vec_double3 vec_make_double3(double x, double y,
                                                     double z)
{
    vec_double3 result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
};
static inline VEC_CFUNC vec_double3 vec_make_double3(vec_double4 other)
{
    vec_double3 result;
    result.xyz = other.xyz;
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
    result.xyz = other;
    result.w = w;
    return result;
}
static inline VEC_CFUNC vec_double4 vec_make_double4(vec_double3 other)
{
    vec_double4 result;
    result.xyz = other;
    return result;
}
