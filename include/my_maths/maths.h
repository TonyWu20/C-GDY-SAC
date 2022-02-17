#pragma once
#if defined __arm64__
#include <arm_neon.h>
#endif
#include "my_maths/base.h"
#include "my_maths/vector.h"

static inline VEC_CFUNC vec_double4 __tg_sqrt(vec_double4 x);
static VEC_CFUNC vec_double2 __tg_sqrt(vec_double2 x)
{
#if defined __SSE2__
    return _mm_sqrt_pd(x);
#elif defined __arm64__
    return vsqrtq_f64(x);
#else
    return vec_make_double2(sqrt(x.x), sqrt(x.y));
#endif
}
static VEC_CFUNC vec_double3 __tg_sqrt(vec_double3 x)
{
    return vec_make_double3(__tg_sqrt(vec_make_double4_undef(x)));
}
static VEC_CFUNC vec_double4 __tg_sqrt(vec_double4 x)
{
#if defined __AVX__
    return _mm256_sqrt_pd(x);
#else
    return vec_make_double4(__tg_sqrt(x.lo), __tg_sqrt(x.hi));
#endif
}
