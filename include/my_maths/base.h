/* This header tries to port simd of macOS
 */
#define VEC_OVERLOAD __attribute__((__overloadable__))
#define VEC_INLINE __attribute__((__always_inline__))
#define VEC_NODEBUG __attribute__((__nodebug__))
#define VEC_CONST __attribute__((__const__))
#define VEC_CFUNC VEC_OVERLOAD VEC_INLINE VEC_NODEBUG VEC_CONST
#define VEC_NOINLINE VEC_OVERLOAD VEC_NODEBUG VEC_CONST
