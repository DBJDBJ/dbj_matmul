#ifndef DBJ_DBJ_MATRIX_INC
#define DBJ_DBJ_MATRIX_INC

#ifndef __clang__
#error please use clang compiler
#error
#else 
#pragma clang system_header 
#endif // __clang__

#include <stdlib.h>
/*
 DBJ MATRIX FRAMEWORK

 (c) 2021 by dbj@dbj.org -- https://dbj.org/license_dbj

usage:

typedef DBJ_MATRIX_STRUCT(char) char_matrix ;

char_matrix * char_matrix_struct_pointer = NULL ;

DBJ_MATRIX_NEW( char_matrix_struct_pointer, rows , cols ) ;

assert(char_matrix_struct_pointer) ;

DBJ_MATRIX_ALIAS(matrix, char, cols) ;

DBJ_MATRIX_CAST(mx, matrix, char_matrix_struct_pointer);

mx[0][0] = '!' ;

*/

#ifndef DBJ_SANITY_MAX_ROW_COUNT
#define DBJ_SANITY_MAX_ROW_COUNT 0xFFFF
#endif

#ifndef DBJ_SANITY_MAX_COL_COUNT
#define DBJ_SANITY_MAX_COL_COUNT 0xFFFF
#endif

#define DBJ_MATRIX_HANDLE void *

#define DBJ_MATRIX_STRUCT(T_) \
struct {\
	unsigned rows;\
	unsigned cols;\
	T_ data[];\
}

#define DBJ_MATRIX_STRUCT_SIZE( T_,R_,C_) \
(sizeof( DBJ_MATRIX_STRUCT(T_) ) + sizeof( T_[R_ * C_] ))

static inline void* DBJ_MATRIX_ALLOC
(const unsigned rows_, const unsigned cols_, size_t type_size_)
{
	#ifndef NDEBUG
	if (rows_ > DBJ_SANITY_MAX_ROW_COUNT) return NULL;
	if (cols_ > DBJ_SANITY_MAX_COL_COUNT) return NULL;
	#endif // NDEBUG

	return calloc(1, type_size_);
}

#define DBJ_MATRIX_NEW(N_,T_,R_,C_) \
do { \
assert(N_ == NULL) ; \
DBJ_MATRIX_STRUCT(T_) * retval = DBJ_MATRIX_ALLOC(R_,C_, DBJ_MATRIX_STRUCT_SIZE(T_,R_,C_) ); \
	 if (retval) { \
		 retval->rows = R_ ; \
		 retval->cols = C_ ; \
	 } ; \
	 N_ = retval ; \
} while(0)

#define DBJ_MATRIX_FREE(M_) do { assert(M_); free(M_); M_ = NULL; } while(0)

#define DBJ_MATRIX_DATA_POINTER(T_, M_) (T_(*)[M_->rows][M_->cols])M_->data

// matrix type is to be used with [][]
#define DBJ_MATRIX_TYPE(T_, C_) T_(*)[C]

#define DBJ_MATRIX_ALIAS(N_,T_,C_) typedef T_(*N_)[C_]

#define DBJ_MATRIX_CAST(N_,A_, M_) A_ N_ = (A_)(M_->data)

#endif // DBJ_DBJ_MATRIX_INC