/*
 Started from systemd macro.h -- SPDX-License-Identifier: LGPL-2.1-or-later

 The rest (c) 2021 July by dbj~dbj.org -- https://dbj.org/license_dbj
*/
#pragma once

#ifndef __clang__
#error Pleae use clang compiler
#endif

#pragma clang system_header

#if __STDC_VERSION__ < 201112L
#error Please use C11 or better
#endif

#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>

#ifndef _WIN32
#include <sys/param.h>
#include <sys/sysmacros.h>
#endif

#include <sys/types.h>

/*
 * STRLEN - return the length of a string literal, minus the trailing NUL byte.
 *          Contrary to strlen(), this is a constant expression.
 * @x: a string literal.
 */
#define STRLEN(x) (sizeof(""x"") - 1)

#define XCONCATENATE(x, y) x##y
#define CONCATENATE(x, y) XCONCATENATE(x, y)

#if defined(static_assert)
#define assert_cc(expr) static_assert(expr, #expr)
#else
#define assert_cc(expr)                                                        \
  struct CONCATENATE(_assert_struct_, __COUNTER__) {                           \
	char x[(expr) ? 0 : -1];                                                   \
  }
#endif

#define _const_ __attribute__((__const__))
#define _pure_ __attribute__((__pure__))
#define _unused_ __attribute__((__unused__))
#define _cleanup_(x) __attribute__((__cleanup__(x)))

 // usage: type_name _auto_cleanup_ pointer_name ;
 // example
 // char * _auto_cleanup_ str = strdup("string content");
 //
#define _auto_cleanup_ _cleanup_(cleanup_free)

 // DBJ 2021-JUL-05
static inline void cleanup_free(void* p) {
	free(*(void**)p);
}

#define _printf_(a, b) __attribute__((__format__(printf, a, b)))
#ifdef __clang__
#  define _alloc_(...)
#else
#  define _alloc_(...) __attribute__((__alloc_size__(__VA_ARGS__)))
#endif
#define _sentinel_ __attribute__((__sentinel__))
#define _section_(x) __attribute__((__section__(x)))
#define _used_ __attribute__((__used__))
#define _destructor_ __attribute__((__destructor__))

// DBJ added 2021 JUL 05
#define _constructor_ __attribute__((__constructor__))

#define _deprecated_ __attribute__((__deprecated__))
#define _packed_ __attribute__((__packed__))
#define _malloc_ __attribute__((__malloc__))
#define _weak_ __attribute__((__weak__))
#define _likely_(x) (__builtin_expect(!!(x), 1))
#define _unlikely_(x) (__builtin_expect(!!(x), 0))
#define _public_ __attribute__((__visibility__("default")))
#define _hidden_ __attribute__((__visibility__("hidden")))
#define _weakref_(x) __attribute__((__weakref__(#x)))
#define _align_(x) __attribute__((__aligned__(x)))
#define _alignas_(x) __attribute__((__aligned__(__alignof(x))))
#define _alignptr_ __attribute__((__aligned__(sizeof(void*))))
#if __GNUC__ >= 7
#define _fallthrough_ __attribute__((__fallthrough__))
#else
#define _fallthrough_
#endif
/* Define C11 noreturn without <stdnoreturn.h> and even on older gcc
 * compiler versions */
#ifndef _noreturn_
#if __STDC_VERSION__ >= 201112L
#define _noreturn_ _Noreturn
#else
#define _noreturn_ __attribute__((__noreturn__))
#endif
#endif

 /* Temporarily disable some warnings */
#define DISABLE_WARNING_DEPRECATED_DECLARATIONS                         \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")

#define DISABLE_WARNING_FORMAT_NONLITERAL                               \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wformat-nonliteral\"")

#define DISABLE_WARNING_MISSING_PROTOTYPES                              \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wmissing-prototypes\"")

#define DISABLE_WARNING_NONNULL                                         \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wnonnull\"")

#define DISABLE_WARNING_SHADOW                                          \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wshadow\"")

#define DISABLE_WARNING_INCOMPATIBLE_POINTER_TYPES                      \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wincompatible-pointer-types\"")

#if HAVE_WSTRINGOP_TRUNCATION
#  define DISABLE_WARNING_STRINGOP_TRUNCATION                           \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wstringop-truncation\"")
#else
#  define DISABLE_WARNING_STRINGOP_TRUNCATION                           \
		_Pragma("GCC diagnostic push")
#endif

#define DISABLE_WARNING_FLOAT_EQUAL \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wfloat-equal\"")

#define DISABLE_WARNING_TYPE_LIMITS \
		_Pragma("GCC diagnostic push");                                 \
		_Pragma("GCC diagnostic ignored \"-Wtype-limits\"")

#define REENABLE_WARNING                                                \
		_Pragma("GCC diagnostic pop")


#define XSTRINGIFY(x) #x
#define STRINGIFY(x) XSTRINGIFY(x)


#define SWAP_TWO(x, y) do {                        \
				typeof(x) _t = (x);                \
				(x) = (y);                         \
				(y) = (_t);                        \
		} while (false)

#define STRV_MAKE(...) ((char**) ((const char*[]) { __VA_ARGS__, NULL }))
#define STRV_MAKE_EMPTY ((char*[1]) { NULL })
#define STRV_MAKE_CONST(...) ((const char* const*) ((const char*[]) { __VA_ARGS__, NULL }))

/* A macro to force copying of a variable from memory. This is useful whenever we want to read something from
 * memory and want to make sure the compiler won't optimize away the destination variable for us. It's not
 * supposed to be a full CPU memory barrier, i.e. CPU is still allowed to reorder the reads, but it is not
 * allowed to remove our local copies of the variables. We want this to work for unaligned memory, hence
 * memcpy() is great for our purposes. */
#define READ_NOW(x)                                                     \
		({                                                              \
				typeof(x) _copy;                                        \
				memcpy(&_copy, &(x), sizeof(_copy));                    \
				asm volatile ("" : : : "memory");                       \
				_copy;                                                  \
		})

 /*
 * ----------------------------------------------------------------------------------------------------
 */
#define DBJ_MATRIX_SIDE_DIMENSION 1024
#define DBJ_MATRIX_DATA_TYPE double
