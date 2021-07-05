
#define DBJ_MATMUL_IMPLEMENTATION
#include "dbj_matmul.h"

#include <time.h>
#include <omp.h>

#define MATMUL_PRINT(...) fprintf( stderr, ## __VA_ARGS__ )

/* Prints vector, or matrix. */
_const_ static
void print(const double* a, const unsigned n_rows_a, const unsigned n_cols_a) {
	for (unsigned i = 0; i < n_rows_a; i++) {
		for (unsigned j = 0; j < n_cols_a; j++) {
			MATMUL_PRINT("%8.3f ", a[i * n_cols_a + j]);
		}
		MATMUL_PRINT("\n");
	}
	MATMUL_PRINT("\n");
}


int main(int argc, char* argv[]) {
	(void)(argc);
	(void)(argv);
	/* Intializes random number generator */
	time_t t;
	srand((unsigned)time(&t));
	srand(0);

	/* For measuring time */
	double t0, t1;

	// const unsigned scale = 10; // provokes matrix print 
	const unsigned scale = 400;
	const unsigned n_rows_a = 4 * scale;
	const unsigned n_cols_a = 3 * scale;
	const unsigned n_rows_b = 3 * scale;
	const unsigned n_cols_b = 2 * scale;

	double* _cleanup_(cleanup_free)  a = NULL;
	double* _cleanup_(cleanup_free)  b = NULL;
	double* _cleanup_(cleanup_free)  c = NULL;
	double* _cleanup_(cleanup_free)  d = NULL;

	ALLOC_WITH_POLICY(a, n_rows_a * n_cols_a * sizeof(*a));
	ALLOC_WITH_POLICY(b, n_rows_b * n_cols_b * sizeof(*b));

	init_rand(a, n_rows_a, n_cols_a);
	init_rand(b, n_rows_b, n_cols_b);

	init_seq(a, n_rows_a, n_cols_a);
	init_seq(b, n_rows_b, n_cols_b);

	t0 = omp_get_wtime();
	c = dot_simple(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b);
	t1 = omp_get_wtime();
	MATMUL_PRINT("Dot Simple: Elapsed time %.3f s\n", t1 - t0);

	t0 = omp_get_wtime();
	d = dot(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b);
	t1 = omp_get_wtime();
	MATMUL_PRINT("Dot: Elapsed time %.3f s\n", t1 - t0);

	if (scale < 11) {
		MATMUL_PRINT("Matrix A:\n");
		print(a, n_rows_a, n_cols_a);
		MATMUL_PRINT("Matrix B:\n");
		print(b, n_rows_b, n_cols_b);
		MATMUL_PRINT("Matrix C:\n");
		print(c, n_rows_a, n_cols_b);
		MATMUL_PRINT("Matrix D:\n");
		print(d, n_rows_a, n_cols_b);
	}

#ifdef WIN32
	system("pause");
#endif
	return EXIT_SUCCESS;
}

