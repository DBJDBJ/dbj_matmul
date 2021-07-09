
#include "ubench.h/ubench.h"

#if 0 // ---------------------------------------------------------------------------------------

#define DBJ_MATMUL_IMPLEMENTATION
#include "dbj_matmul.h"

#define DBJ_MATMUL_OMP_IMP
#include "dbj_matmul_omp.h"

#undef DBJ_API
#define DBJ_API __attribute__((const)) static

//////////////////////////////////////////////////////////////////
// ubench bench functions have no parameters
// thus we use common data aka globals

static struct app_data_type {

	const unsigned n_rows_a;
	const unsigned n_cols_a;
	const unsigned n_rows_b;
	const unsigned n_cols_b;

	unsigned start_counter;

	// the matrixes
	double* a;
	double* b;
	double* rezult; /* rezult size is a rows * b cols */

} app_data = {
	.n_rows_a = DBJ_MATRIX_SIDE_DIMENSION ,
	.n_cols_a = DBJ_MATRIX_SIDE_DIMENSION ,
	.n_rows_b = DBJ_MATRIX_SIDE_DIMENSION ,
	.n_cols_b = DBJ_MATRIX_SIDE_DIMENSION
	// the rest is auto zeroed in C
};

DBJ_API void app_start(void) _constructor_;
DBJ_API void app_end(void) _destructor_;

DBJ_API void app_start(void)
{
	ALLOC_WITH_POLICY(app_data.a, app_data.n_rows_a * app_data.n_cols_a * sizeof(*app_data.a));
	ALLOC_WITH_POLICY(app_data.b, app_data.n_rows_b * app_data.n_cols_b * sizeof(*app_data.b));
	ALLOC_WITH_POLICY(app_data.rezult, app_data.n_rows_a * app_data.n_cols_b * sizeof(*app_data.rezult));

	init_seq(app_data.a, app_data.n_rows_a, app_data.n_cols_a);
	init_seq(app_data.b, app_data.n_rows_b, app_data.n_cols_b);
	init_seq(app_data.rezult, app_data.n_rows_a, app_data.n_cols_b);
}

DBJ_API void app_end(void)
{
	free(app_data.a);
	free(app_data.b);
	free(app_data.rezult);
}

UBENCH(dbj_, simple) {
	dot_simple(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b, app_data.rezult);
}

UBENCH(dbj_, faster) {
	dot(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b, app_data.rezult);
}

UBENCH(dbj_, omp_simple) {
	omp_dot_simple(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b, app_data.rezult);
}

UBENCH(dbj_, omp_faster) {
	omp_dot_faster(app_data.a, app_data.n_rows_a, app_data.n_cols_a, app_data.b, app_data.n_rows_b, app_data.n_cols_b, app_data.rezult);
}
#endif // 0 ---------------------------------------------------------------------------------------

UBENCH_STATE();

int main(const int argc, const char** argv)
{
	(void)argc;
	(void)argv;
	// WIN32 hack to kick-start the ANSI ESC codes interpreter
#ifdef WIN32
	system(" ");
#endif
	return ubench_main(argc, argv);
}


