#pragma once
#include "my_maths/base.h"
#include "my_maths/vector.h"

static vec_double3x3 VEC_CFUNC vec_diagonal_matrix(vec_double3 __x);
static vec_double4x4 VEC_CFUNC vec_diagonal_matrix(vec_double4 __x);
#define matrix_from_diagonal vec_diagonal_matrix

// clang-format off
static vec_double3x3 VEC_CFUNC vec_diagonal_matrix(vec_double3 __x) { vec_double3x3 __r = { .columns[0] = {__x.x,0,0}, .columns[1] = {0,__x.y,0}, .columns[2] = {0,0,__x.z} }; return __r; }
static vec_double4x4 VEC_CFUNC vec_diagonal_matrix(vec_double4 __x) { vec_double4x4 __r = { .columns[0] = {__x.x,0,0,0}, .columns[1] = {0,__x.y,0,0}, .columns[2] = {0,0,__x.z,0}, .columns[3] = {0,0,0,__x.w} }; return __r; }
// clang-format on

static vec_double3x3 VEC_CFUNC vec_matrix(vec_double3 col0, vec_double3 col1,
                                          vec_double3 col2);
static vec_double4x4 VEC_CFUNC vec_matrix(vec_double4 col0, vec_double4 col1,
                                          vec_double4 col2, vec_double4 col3);
#define matrix_from_columns vec_matrix

// clang-format off
static vec_double3x3 VEC_CFUNC vec_matrix(vec_double3 col0, vec_double3 col1, vec_double3 col2) { vec_double3x3 __r = { .columns[0] = col0, .columns[1] = col1, .columns[2] = col2 }; return __r; }
static vec_double4x4 VEC_CFUNC vec_matrix(vec_double4 col0, vec_double4 col1, vec_double4 col2, vec_double4 col3) { vec_double4x4 __r = { .columns[0] = col0, .columns[1] = col1, .columns[2] = col2, .columns[3] = col3 }; return __r; }
// clang-format on

static vec_double3x3 VEC_NOINLINE vec_matrix3x3(vec_quatd q);
static vec_double4x4 VEC_NOINLINE vec_matrix4x4(vec_quatd q);

static vec_double3x3 VEC_NOINLINE vec_matrix3x3(vec_quatd q)
{
    vec_double4x4 r = vec_matrix4x4(q);
    return (vec_double3x3){r.columns[0].xyz, r.columns[1].xyz,
                           r.columns[2].xyz};
}

static vec_double4x4 VEC_NOINLINE vec_matrix4x4(vec_quatd q)
{
    vec_double4 v = q.vector;
    vec_double4x4 r = {.columns[0] = {1 - 2 * (v.y * v.y + v.z * v.z),
                                      2 * (v.x * v.y + v.z * v.w),
                                      2 * (v.x * v.z - v.y * v.w), 0},
                       .columns[1] = {2 * (v.x * v.y - v.z * v.w),
                                      1 - 2 * (v.z * v.z + v.x * v.x),
                                      2 * (v.y * v.z + v.x * v.w), 0},
                       .columns[2] = {2 * (v.z * v.x + v.y * v.w),
                                      2 * (v.y * v.z - v.x * v.w),
                                      1 - 2 * (v.y * v.y + v.x * v.x), 0},
                       .columns[3] = {0, 0, 0, 1}};
    return r;
}
// clang-format off
static vec_double3x3 VEC_CFUNC vec_mul(double __a, vec_double3x3 __x) { __x.columns[0] *= __a; __x.columns[1] *= __a; __x.columns[2] *= __a; return __x; }
static vec_double4x4 VEC_CFUNC vec_mul(double __a, vec_double4x4 __x) { __x.columns[0] *= __a; __x.columns[1] *= __a; __x.columns[2] *= __a; __x.columns[3] *= __a; return __x; }
static vec_double3 VEC_CFUNC vec_mul(vec_double3x3 __x, vec_double3 __y) { vec_double3 __r = __x.columns[0]*__y[0]; __r = vec_muladd( __x.columns[1], __y[1],__r); __r = vec_muladd( __x.columns[2], __y[2],__r); return __r; }
static vec_double4 VEC_CFUNC vec_mul(vec_double4x4 __x, vec_double4 __y) { vec_double4 __r = __x.columns[0]*__y[0]; __r = vec_muladd( __x.columns[1], __y[1],__r); __r = vec_muladd( __x.columns[2], __y[2],__r); __r = vec_muladd( __x.columns[3], __y[3],__r); return __r; }
static vec_double3 VEC_CFUNC vec_mul(vec_double3 __x, vec_double3x3 __y) { return vec_mul(vec_transpose(__y), __x); }
static vec_double4 VEC_CFUNC vec_mul(vec_double4 __x, vec_double4x4 __y) { return vec_mul(vec_transpose(__y), __x); }
// clang-format on
