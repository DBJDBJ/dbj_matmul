
#define DBJ_MATMUL_IMPLEMENTATION
#include "dbj_matmul.h"

#define DBJ_MATMUL_OMP_IMP
#include "dbj_matmul_omp.h"

#include "ubench.h/ubench.h"

//////////////////////////////////////////////////////////////////
// ubench bench functions have no parameters
// thus we use common data aka globals

#define SCALE 10U
// #define SCALE 400U

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

		system(" "); // WIN32 hack to kick-start the ANSI ESC codes interpreter
	}

	app_data.start_counter += 1;
}

static void app_end(void)
{
	free(app_data.a);
	free(app_data.b);
}

UBENCH(dbj_matmul, dot_simple) {
	app_data.c = dot_simple(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	free(app_data.c);
}

UBENCH(dbj_matmul, dot_faster) {
	app_data.d = dot(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	free(app_data.d);
}

UBENCH(dbj_matmul, omp_dot_simple) {
	app_data.c = omp_dot_simple(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	free(app_data.c);
}

UBENCH(dbj_matmul, omp_dot_faster) {
	app_data.d = omp_dot_faster(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	free(app_data.d);
}

UBENCH_MAIN();


