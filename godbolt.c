/////////////////////////////////////////////////////////////////////////
// https://godbolt.org/z/7nqTx7nxY
/////////////////////////////////////////////////////////////////////////

#include "dbj_matmul_common.h"
#define DBJ_MATMUL_IMPLEMENTATION
#include "dbj_matmul.h"
#include "ubench.h/ubench.h"


/////////////////////////////////////////////////////////////////////////
#define DBJ_MATMUL_IMPLEMENTATION
#ifdef DBJ_MATMUL_IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
typedef double dbj_matrix_data_type;
/////////////////////////////////////////////////////////////////////////

#if 0 // --------------------------------------------------------------------------------


/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
	It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
	of m, meaning, rows_t == cols_m, and cols_t == rows_m.
	The original matrix m stays intact. */
DBJ_API
dbj_matrix_data_type* transpose(
	const dbj_matrix_data_type* m,
	const unsigned rows_m, const unsigned cols_m,
	dbj_matrix_data_type* t) {
	for (size_t i = 0; i < rows_m; i++) {
		for (size_t j = 0; j < cols_m; j++) {
			t[j * rows_m + i] = m[i * cols_m + j];
		}
	}

	return t;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant doesn't transpose matrix b, and it's a lot slower. */
DBJ_API
dbj_matrix_data_type* dot_simple(const dbj_matrix_data_type* a, const unsigned rows_a, const unsigned cols_a,
	const dbj_matrix_data_type* b, const unsigned rows_b, const unsigned cols_b, dbj_matrix_data_type* c) {

	assert(cols_a == rows_b);

	for (size_t i = 0; i < rows_a; i++) {
		for (size_t k = 0; k < cols_b; k++) {
			dbj_matrix_data_type sum = 0.0;
			for (size_t j = 0; j < cols_a; j++) {
				sum += a[i * cols_a + j] * b[j * cols_b + k];
			}
			c[i * cols_b + k] = sum;
		}
	}

	return c;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant transposes matrix b, and it's a lot faster. */
DBJ_API
dbj_matrix_data_type* dot(const dbj_matrix_data_type* a, const unsigned rows_a, const unsigned cols_a,
	const dbj_matrix_data_type* b, const unsigned rows_b, const unsigned cols_b, dbj_matrix_data_type* c) {

	assert(cols_a == rows_b);

	dbj_matrix_data_type* bt = 0;
	ALLOC_WITH_POLICY(bt, rows_b * cols_b * sizeof(*b));

	bt = transpose(b, rows_b, cols_b, bt);

	for (unsigned i = 0; i < rows_a; i++) {
		for (unsigned k = 0; k < cols_b; k++) {
			dbj_matrix_data_type sum = 0.0;
			for (unsigned j = 0; j < cols_a; j++) {
				sum += a[i * cols_a + j] * bt[k * rows_b + j];
			}
			c[i * cols_b + k] = sum;
		}
	}
	return c;
}

DBJ_API
void init_seq(dbj_matrix_data_type* a, const unsigned rows_a, const unsigned cols_a) {
	for (unsigned i = 0; i < rows_a; i++) {
		for (unsigned j = 0; j < cols_a; j++) {
			a[i * cols_a + j] = i * cols_a + j;
		}
	}
}

#endif // 0 --------------------------------------------------------------------------------

DBJ_API void* matmul_surprise(
	const unsigned a_rows,
	const unsigned cols_a,
	const unsigned b_rows,
	const unsigned b_cols,
	const unsigned c_rows,
	const unsigned c_cols,
	dbj_matrix_data_type(*ax)[cols_a],
	dbj_matrix_data_type(*bx)[b_cols],
	dbj_matrix_data_type(*mx)[c_cols])
{
	assert(b_rows == cols_a);
	assert(a_rows == c_rows);
	assert(b_cols == c_cols);

	// VMT is runtime mechanism
	// runtime casting takes time
	// matridbj_x_data_type *ax)[a_cols] = a;
	// matridbj_x_data_type *bx)[b_cols] = b;
	// matridbj_x_data_type *mx)[b_rows] = m;

	for (unsigned i = 0; i < a_rows; ++i)
	{
		for (unsigned j = 0; j < b_cols; ++j)
		{
			dbj_matrix_data_type t = 0.0;
			for (unsigned k = 0; k < cols_a; ++k)
				t += ax[i][k] * bx[k][j];
			mx[i][j] = t;
		}
	}
	return mx;
}

DBJ_API void* simple_mat_mul_0(
	const unsigned a_rows, const unsigned a_cols,
	const unsigned b_rows, const unsigned b_cols,
	const unsigned m_rows, const unsigned m_cols,
	dbj_matrix_data_type a[static a_rows][a_cols],
	dbj_matrix_data_type b[static b_rows][b_cols],
	dbj_matrix_data_type m[static m_rows][m_cols])
{
	assert(b_rows == a_cols);
	assert(a_rows == m_rows);
	assert(b_cols == m_cols);

	for (unsigned i = 0; i < a_rows; ++i)
	{
		for (unsigned j = 0; j < b_cols; ++j)
		{
			dbj_matrix_data_type t = 0.0;
			for (unsigned k = 0; k < a_cols; ++k)
				t += a[i][k] * b[k][j];
			m[i][j] = t;
		}
	}
	return m;
}
/////////////////////////////////////////////////////////////////////////
#endif // DBJ_MATMUL_IMPLEMENTATION

//////////////////////////////////////////////////////////////////
// ubench bench functions have no parameters
// thus we use common data aka globals

#define DBJ_MATRIX_SIDE_DIMENSION 0xF

static struct app_data_type {

	const unsigned rows_a;
	const unsigned cols_a;
	const unsigned rows_b;
	const unsigned cols_b;
	const unsigned rows_r;
	const unsigned cols_r;
	// the matrixes
	dbj_matrix_data_type* a;
	dbj_matrix_data_type* b;
	dbj_matrix_data_type* r; /* rezult size is a rows * b cols */

} app_data = {
	.rows_a = DBJ_MATRIX_SIDE_DIMENSION ,
	.cols_a = DBJ_MATRIX_SIDE_DIMENSION ,
	.rows_b = DBJ_MATRIX_SIDE_DIMENSION ,
	.cols_b = DBJ_MATRIX_SIDE_DIMENSION ,
	.rows_r = DBJ_MATRIX_SIDE_DIMENSION , /* A rows */
	.cols_r = DBJ_MATRIX_SIDE_DIMENSION   /* B cols */
	// the rest is auto zeroed 
};

__attribute__((constructor)) DBJ_API void app_start(void)
{
	assert(app_data.rows_b == app_data.cols_a);
	assert(app_data.rows_a == app_data.rows_r);
	assert(app_data.cols_b == app_data.cols_r);

	ALLOC_WITH_POLICY(app_data.a, app_data.rows_a * app_data.cols_a * sizeof(*app_data.a));
	ALLOC_WITH_POLICY(app_data.b, app_data.rows_b * app_data.cols_b * sizeof(*app_data.b));
	ALLOC_WITH_POLICY(app_data.r, app_data.rows_a * app_data.cols_b * sizeof(*app_data.r));

	init_seq(app_data.a, app_data.rows_a, app_data.cols_a);
	init_seq(app_data.b, app_data.rows_b, app_data.cols_b);
	init_seq(app_data.r, app_data.rows_a, app_data.cols_b);
}

__attribute__((destructor)) DBJ_API void app_end(void)
{
	free(app_data.a);
	free(app_data.b);
	free(app_data.r);
}
/////////////////////////////////////////////////////////////////////////

UBENCH(godbolt, simple) {
	dot_simple(app_data.a, app_data.rows_a, app_data.cols_a, app_data.b, app_data.rows_b, app_data.cols_b, app_data.r);
}

UBENCH(godbolt, not_faster) {
	dot(app_data.a, app_data.rows_a, app_data.cols_a, app_data.b, app_data.rows_b, app_data.cols_b, app_data.r);
}

UBENCH(godbolt, snazzy) {
	matmul_surprise(
		app_data.rows_a,
		app_data.cols_a,
		app_data.rows_b,
		app_data.cols_b,
		app_data.rows_r,
		app_data.cols_r,
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r);
}

UBENCH(godbolt, simple_mat_mul_0) {
	simple_mat_mul_0(
		app_data.rows_a,
		app_data.cols_a,
		app_data.rows_b,
		app_data.cols_b,
		app_data.rows_r,
		app_data.cols_r,
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r);
}

/////////////////////////////////////////////////////////////////////////


