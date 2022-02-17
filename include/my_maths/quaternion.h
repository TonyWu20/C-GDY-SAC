#pragma once
#include "my_maths/base.h"
#include "my_maths/geometry.h"
#include "my_maths/vector.h"
#include <math.h>

static inline VEC_CFUNC vec_double3 vec_imag(vec_quatd q)
{
    return q.vector.xyz;
}
static inline VEC_CFUNC double vec_real(vec_quatd q)
{
    return q.vector.w;
}
static inline VEC_CFUNC vec_quatd vec_quaternion(double ix, double iy,
                                                 double iz, double r)
{
    return (vec_quatd){{ix, iy, iz, r}};
}
static inline VEC_CFUNC vec_quatd vec_quaternion(vec_double4 xyzr)
{
    return (vec_quatd){xyzr};
}
static inline VEC_CFUNC vec_quatd _vec_quaternion(vec_double3 imag, double real)
{
    return vec_quaternion(vec_make_double4(imag, real));
}
static inline VEC_CFUNC vec_quatd _vec_quaternion_reduced(vec_double3 from,
                                                          vec_double3 to)
{
    vec_double3 half = vec_normalize(from + to);
    return _vec_quaternion(vec_cross(from, half), vec_dot(from, half));
}
static inline VEC_CFUNC vec_quatd vec_quaternion(double angle,
                                                 vec_double3 axis);
static inline VEC_CFUNC vec_double3 vec_act(vec_quatd q, vec_double3 v);
static inline VEC_CFUNC vec_quatd vec_quaternion(double angle, vec_double3 axis)
{
    return _vec_quaternion(sin(angle / 2) * axis, cos(angle / 2));
}
static inline VEC_CFUNC vec_double3 vec_act(vec_quatd q, vec_double3 v)
{
    vec_double3 t = 2 * vec_cross(vec_imag(q), v);
    return v + vec_real(q) * t + vec_cross(vec_imag(q), t);
}
