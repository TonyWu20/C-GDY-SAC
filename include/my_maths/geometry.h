#pragma once
#include "my_maths/base.h"
#include "my_maths/common.h"
#include "my_maths/vector.h"

// clang-format off
static double VEC_CFUNC vec_dot(vec_double3 __x, vec_double3 __y);
static double VEC_CFUNC vec_dot(vec_double4 __x, vec_double4 __y);
#define vector_dot vec_dot
// clang-format on
//
static double VEC_CFUNC vec_length(vec_double3 __x);
static double VEC_CFUNC vec_length(vec_double4 __x);
#define vector_length vec_length

static double VEC_CFUNC vec_length_squared(vec_double3 __x);
static double VEC_CFUNC vec_length_squared(vec_double4 __x);
#define vector_length_squared vec_length_squared

static vec_double3 VEC_CFUNC vec_normalize(vec_double3 __x);
static vec_double4 VEC_CFUNC vec_normalize(vec_double4 __x);
#define vector_normalize vec_normalize

static vec_double3 VEC_CFUNC vec_cross(vec_double3 __x, vec_double3 __y);
#define vector_cross vec_cross

#pragma mark - Implementation

static double VEC_CFUNC vec_dot(vec_double3 __x, vec_double3 __y)
{
    return vec_reduce_add(__x * __y);
}
static double VEC_CFUNC vec_dot(vec_double4 __x, vec_double4 __y)
{
    return vec_reduce_add(__x * __y);
}

static double VEC_CFUNC vec_precise_length(vec_double3 __x)
{
    return sqrt(vec_length_squared(__x));
}
static double VEC_CFUNC vec_precise_length(vec_double4 __x)
{
    return sqrt(vec_length_squared(__x));
}
static double VEC_CFUNC vec_length_squared(vec_double3 __x)
{
    return vec_dot(__x, __x);
}
static double VEC_CFUNC vec_length_squared(vec_double4 __x)
{
    return vec_dot(__x, __x);
}
static double VEC_CFUNC vec_length(vec_double3 __x)
{
    return vec_precise_length(__x);
}
static double VEC_CFUNC vec_length(vec_double4 __x)
{
    return vec_precise_length(__x);
}

static vec_double3 VEC_CFUNC vec_precise_normalize(vec_double3 __x)
{
    return __x * vec_precise_rsqrt(vec_length_squared(__x));
}
static vec_double4 VEC_CFUNC vec_precise_normalize(vec_double4 __x)
{
    return __x * vec_precise_rsqrt(vec_length_squared(__x));
}
static vec_double3 VEC_CFUNC vec_normalize(vec_double3 __x)
{
    return vec_precise_normalize(__x);
}
static vec_double4 VEC_CFUNC vec_normalize(vec_double4 __x)
{
    return vec_precise_normalize(__x);
}

static vec_double3 VEC_CFUNC vec_cross(vec_double3 __x, vec_double3 __y)

{
    return (__x.zxy * __y - __x * __y.zxy).zxy;
}
