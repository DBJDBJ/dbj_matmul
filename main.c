
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

//////////////////////////////////////////////////////////////////
/// ubench function have no parameters
/// thus we use common data aka globals

#define SCALE 400U

static struct app_data_type {

	const unsigned scale;
	const unsigned n_rows_a; // 4 * scale;
	const unsigned n_cols_a; // 3 * scale;
	const unsigned n_rows_b; // 3 * scale;
	const unsigned n_cols_b; // 2 * scale;

	unsigned start_counter;

	// the matrixes
	double* a;
	double* b;
	double* c;
	double* d;

} app_data = {
	.scale = SCALE ,
	.n_rows_a = 4U * SCALE ,
	.n_cols_a = 3U * SCALE ,
	.n_rows_b = 3U * SCALE ,
	.n_cols_b = 2U * SCALE
	// the rest is auto zeroed in C
};

static void app_start(void) _constructor_;
static void app_end(void) _destructor_;

static void app_start(void)
{
	if (app_data.start_counter < 1) {
		ALLOC_WITH_POLICY(app_data.a, app_data.n_rows_a * app_data.n_cols_a * sizeof(*app_data.a));
		ALLOC_WITH_POLICY(app_data.b, app_data.n_rows_b * app_data.n_cols_b * sizeof(*app_data.b));

		init_rand(app_data.a, app_data.n_rows_a, app_data.n_cols_a);
		init_rand(app_data.b, app_data.n_rows_b, app_data.n_cols_b);

		init_seq(app_data.a, app_data.n_rows_a, app_data.n_cols_a);
		init_seq(app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	}

	app_data.start_counter += 1;
}

static void app_end(void)
{
	free(app_data.a);
	free(app_data.b);
	free(app_data.c);
	free(app_data.d);
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

	t0 = omp_get_wtime();
	app_data.c = dot_simple(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	t1 = omp_get_wtime();
	MATMUL_PRINT("Dot Simple: Elapsed time %.3f s\n", t1 - t0);

	t0 = omp_get_wtime();
	app_data.d = dot(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	t1 = omp_get_wtime();
	MATMUL_PRINT("Dot: Elapsed time %.3f s\n", t1 - t0);

#ifdef WIN32
	system("pause");
#endif
	return EXIT_SUCCESS;
}

