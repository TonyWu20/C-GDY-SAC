#pragma once
#include "my_maths/base.h"
#include "my_maths/maths.h"
#include "my_maths/vector.h"
#include <math.h>
static inline VEC_CFUNC double vec_reduce_add(vec_double2 x);
static inline VEC_CFUNC double vec_reduce_add(vec_double3 x);
static inline VEC_CFUNC double vec_reduce_add(vec_double4 x);

static inline VEC_CFUNC double vec_reduce_add(vec_double2 x)
{
    return x.x + x.y;
}
static inline VEC_CFUNC double vec_reduce_add(vec_double3 x)
{
    return x.x + x.y + x.z;
}
static inline VEC_CFUNC double vec_reduce_add(vec_double4 x)
{
    return vec_reduce_add(x.lo + x.hi);
}

static inline VEC_CFUNC vec_double3 vec_precise_rsqrt(vec_double3 x);
static inline VEC_CFUNC vec_double4 vec_precise_rsqrt(vec_double4 x);

static inline VEC_CFUNC double vec_precise_rsqrt(double x)
{
    return 1 / sqrt(x);
}
static inline VEC_CFUNC vec_double3 vec_precise_rsqrt(vec_double3 x)
{
    return 1 / __tg_sqrt(x);
}

static inline VEC_CFUNC vec_double4 vec_precise_rsqrt(vec_double4 x)
{
    return 1 / __tg_sqrt(x);
}
