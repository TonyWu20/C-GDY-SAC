#pragma once
#include <math.h>
#define PI (atan(1) * 4)
#define VEC_OVERLOAD __attribute__((__overloadable__))
#define VEC_INLINE __attribute__((__always_inline__))
#define VEC_NODEBUG __attribute__((__nodebug__))
#define VEC_CONST __attribute__((__const__))
#define VEC_CFUNC VEC_OVERLOAD VEC_INLINE VEC_NODEBUG VEC_CONST
#define VEC_NOINLINE VEC_OVERLOAD VEC_NODEBUG VEC_CONST
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
