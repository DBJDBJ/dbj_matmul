#ifndef DBJ_MATMUL_INC
#define DBJ_MATMUL_INC


#include "dbj_matmul_common.h"

/*
you know the drill, define DBJ_MATMUL_IMPLEMENTATION
on exactly one C/C++, not two or more
*/


// works server side too
#define NOMEM_POLICY( BOOLEXP_ )\
	if (! BOOLEXP_ ) {\
		perror( __FILE__ ", Could not allocate memory!");\
		exit(-1);\
	}

#define ALLOC_WITH_POLICY(PTR_ , SIZE_)    \
do { \
PTR_ = calloc(1,SIZE_);\
NOMEM_POLICY(PTR_) ;\
} while(0)

/* Initializes vector or matrix, sequentially, with indices. */
_const_
void init_seq(double*, const unsigned, const unsigned);
/* Initializes vector or matrix, randomly. */
_const_
void init_rand(double*, const unsigned, const unsigned);
/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
	It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
	of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
	The original matrix m stays intact. */
_const_
double* transpose(const double*, const unsigned, const unsigned, double*);

/////////////////////////////////////////////////////////////////////////

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant doesn't transpose matrix b, and it's a lot slower. */
_const_
double* dot_simple(const double*, const unsigned, const unsigned,
	const double*, const unsigned, const unsigned, double* /* rezult */);

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant transposes matrix b, and it's a lot faster. */
_const_
double* dot(const double* a, const unsigned n_rows_a, const unsigned n_cols_a, \
	const double* b, const unsigned n_rows_b, const unsigned n_cols_b, double* /* rezult */);

/* Prints vector, or matrix. */
_const_
void print(const double* a, const unsigned n_rows_a, const unsigned n_cols_a);
/////////////////////////////////////////////////////////////////////////
#ifdef DBJ_MATMUL_IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////
/* Initializes vector or matrix, sequentially, with values of indices. */
_const_
void init_seq(double* a, const unsigned n_rows_a, const unsigned n_cols_a) {
	for (unsigned i = 0; i < n_rows_a; i++) {
		for (unsigned j = 0; j < n_cols_a; j++) {
			a[i * n_cols_a + j] = i * n_cols_a + j;
		}
	}
}

/* Initializes vector or matrix, randomly. */
_const_
void init_rand(double* a, const unsigned n_rows_a, const unsigned n_cols_a) {
	for (size_t i = 0; i < n_rows_a; i++) {
		for (size_t j = 0; j < n_cols_a; j++) {
			a[i * n_cols_a + j] = rand() / (double)RAND_MAX;
		}
	}
}

/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
	It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
	of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
	The original matrix m stays intact. */
_const_
double* transpose(const double* m, const unsigned n_rows_m, const unsigned n_cols_m, double* t) {
	for (size_t i = 0; i < n_rows_m; i++) {
		for (size_t j = 0; j < n_cols_m; j++) {
			t[j * n_rows_m + i] = m[i * n_cols_m + j];
		}
	}

	return t;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant doesn't transpose matrix b, and it's a lot slower. */
_const_
double* dot_simple(const double* a, const unsigned n_rows_a, const unsigned n_cols_a,
	const double* b, const unsigned n_rows_b, const unsigned n_cols_b, double* c) {

	assert(n_cols_a == n_rows_b);

	for (size_t i = 0; i < n_rows_a; i++) {
		for (size_t k = 0; k < n_cols_b; k++) {
			double sum = 0.0;
			for (size_t j = 0; j < n_cols_a; j++) {
				sum += a[i * n_cols_a + j] * b[j * n_cols_b + k];
			}
			c[i * n_cols_b + k] = sum;
		}
	}

	return c;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant transposes matrix b, and it's a lot faster. */
_const_
double* dot(const double* a, const unsigned n_rows_a, const unsigned n_cols_a,
	const double* b, const unsigned n_rows_b, const unsigned n_cols_b, double* c) {

	assert(n_cols_a == n_rows_b);

	double* bt = 0;
	ALLOC_WITH_POLICY(bt, n_rows_b * n_cols_b * sizeof(*b));

	bt = transpose(b, n_rows_b, n_cols_b, bt);

	for (unsigned i = 0; i < n_rows_a; i++) {
		for (unsigned k = 0; k < n_cols_b; k++) {
			double sum = 0.0;
			for (unsigned j = 0; j < n_cols_a; j++) {
				sum += a[i * n_cols_a + j] * bt[k * n_rows_b + j];
			}
			c[i * n_cols_b + k] = sum;
		}
	}
	return c;
}

/////////////////////////////////////////////////////////////////////////
#endif // DBJ_MATMUL_IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////

#endif // DBJ_MATMUL_INC
