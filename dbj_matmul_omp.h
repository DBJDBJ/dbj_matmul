#ifndef DBJ_MATMUL_OPEN_MP
#define DBJ_MATMUL_OPEN_MP

/*

Inspiration: Source: https://github.com/ivanbgd/Matrix-Multiplication-MatMul-C
The rest: (c) 2021 by dbj at dbj dot org -- https://dbj.org/license_dbj

Matrices are represented as 1-D arrays in memory.
That means they are contiguous in memory.
Minimum dimension is 1, not 0, and internal dimensions must match.

*/

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/////////////////////////////////////////////////////////////////////////////////

void omp_init_seq(double* a, const unsigned n_rows_a, const unsigned n_cols_a);

void omp_init_rand(double* a, const unsigned n_rows_a, const unsigned n_cols_a);

double* omp_transpose(const double* m, const unsigned n_rows_m, const unsigned n_cols_m, double* t);

double* omp_dot_simple(const double* a, const unsigned n_rows_a, const unsigned n_cols_a,
	const double* b, const unsigned n_rows_b, const unsigned n_cols_b);

double* omp_dot_faster(const double* a, const unsigned n_rows_a, const unsigned n_cols_a,
	const double* b, const unsigned n_rows_b, const unsigned n_cols_b);
/////////////////////////////////////////////////////////////////////////////////

#ifdef DBJ_MATMUL_OMP_IMP
/* Initializes vector or matrix, sequentially, with indices. */
void omp_init_seq(double* a, const unsigned n_rows_a, const unsigned n_cols_a) {
	int i;

#pragma omp parallel for default(none) private(i) shared(a, n_rows_a, n_cols_a) schedule(static)
	for (i = 0; i < n_rows_a; i++) {
		for (size_t j = 0; j < n_cols_a; j++) {
			a[i * n_cols_a + j] = i * n_cols_a + j;
		}
	}
}

/* Initializes vector or matrix, randomly. */
void omp_init_rand(double* a, const unsigned n_rows_a, const unsigned n_cols_a) {
	int i;

	/* Schedule has to be either guided or dynamic; if it's static or runtime, the random numbers repeat. */
#pragma omp parallel for default(none) private(i) shared(a, n_rows_a, n_cols_a) schedule(dynamic)
	for (i = 0; i < n_rows_a; i++) {
		for (size_t j = 0; j < n_cols_a; j++) {
			a[i * n_cols_a + j] = rand() / (double)RAND_MAX;
		}
	}
}

/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
	It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
	of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
	The original matrix m stays intact. */
double* omp_transpose(const double* m, const unsigned n_rows_m, const unsigned n_cols_m, double* t) {
	int i, j;

#pragma omp parallel for default(none) private(i, j) shared(m, n_rows_m, n_cols_m, t) schedule(static)
	for (i = 0; i < n_rows_m; i++) {
		for (j = 0; j < n_cols_m; j++) {
			t[j * n_rows_m + i] = m[i * n_cols_m + j];
		}
	}

	return t;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant doesn't transpose matrix b, and it's a lot slower. */
double* omp_dot_simple(const double* a, const unsigned n_rows_a, const unsigned n_cols_a,
	const double* b, const unsigned n_rows_b, const unsigned n_cols_b) {

	if (n_cols_a != n_rows_b) {
		printf("#columns A must be equal to #rows B!\n");
		system("pause");
		exit(-2);
	}

	double* c = malloc(n_rows_a * n_cols_b * sizeof(*c));
	if (c == NULL) {
		printf("Couldn't allocate memory!\n");
		system("pause");
		exit(-1);
	}

	int i, j, k;

#pragma omp parallel for default(none) private(i, j, k) shared(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c) schedule(static)
	for (i = 0; i < n_rows_a; i++) {
		for (k = 0; k < n_cols_b; k++) {
			double sum = 0.0;
			for (j = 0; j < n_cols_a; j++) {
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
double* omp_dot_faster(const double* a, const unsigned n_rows_a, const unsigned n_cols_a,
	const double* b, const unsigned n_rows_b, const unsigned n_cols_b) {

	int i, j, k;

	if (n_cols_a != n_rows_b) {
		printf("#columns A must be equal to #rows B!\n");
		system("pause");
		exit(-2);
	}

	double* bt = malloc(n_rows_b * n_cols_b * sizeof(*b));

	double* c = malloc(n_rows_a * n_cols_b * sizeof(*c));

	if ((c == NULL) || (bt == NULL)) {
		printf("Couldn't allocate memory!\n");
		system("pause");
		exit(-1);
	}

	bt = omp_transpose(b, n_rows_b, n_cols_b, bt);

#pragma omp parallel for default(none) private(i, j, k) shared(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b, c, bt) schedule(static)
	for (i = 0; i < n_rows_a; i++) {
		for (k = 0; k < n_cols_b; k++) {
			double sum = 0.0;
			for (j = 0; j < n_cols_a; j++) {
				sum += a[i * n_cols_a + j] * bt[k * n_rows_b + j];
			}
			c[i * n_cols_b + k] = sum;
		}
	}

	free(bt);

	return c;
}

#endif // DBJ_MATMUL_OMP_IMP

#endif // DBJ_MATMUL_OPEN_MP
