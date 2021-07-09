/*

 This is benchmarking of a collection of matrix multiplication algorithms.
 Algorithms are kept as simple as possible. No structs are passed as arguments.

 Different comoilers multiplied with different platforms multiplied selection
 of data types  yield a complex picture of
 benchmarking results.

 Use this file to recompile and re measure whenver selecting
 the right matrix multiplication algorithm

 https://godbolt.org/z/joMo4Tn3T

 */



 // #include "https://raw.githubusercontent.com/sheredom/ubench.h/master/ubench.h"
#include "ubench.h/ubench.h"

#ifdef __linux__
#include <x86intrin.h> // #define __SSE__ 1 for SEE versions
#else
#include <intrin.h> // #define __SSE__ 1 for SEE versions
#endif

/* NDEBUG == RELEASE */
#include <assert.h>
/////////////////////////////////////////////////////////////////////////
#define DBJ_MATMUL_IMPLEMENTATION
#ifdef DBJ_MATMUL_IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////

#undef NOMEM_POLICY

#ifdef NDEBUG 
#define NOMEM_POLICY( BOOLEXP_ ) ((void)BOOLEXP_)
#else // ! NDEBUG == DEBUG
#define NOMEM_POLICY( BOOLEXP_ ) if (! BOOLEXP_ ) { perror( __FILE__ ", Could not allocate memory!"); exit(-1); }
#endif // ! NDEBUG

#undef ALLOC_WITH_POLICY
#define ALLOC_WITH_POLICY(PTR_ , SIZE_)  do { PTR_ = calloc(1,SIZE_); NOMEM_POLICY(PTR_); } while(0)

#undef DBJ_API
#define DBJ_API static
/////////////////////////////////////////////////////////////////////////

typedef double dbj_matrix_data_type;
#define dbj_matrix_data_type_name "double"

#define DBJ_MATRIX_SIDE_DIMENSION 0xFF

/////////////////////////////////////////////////////////////////////////
#if __SSE__

DBJ_API void* simple_mat_transpose(
	const unsigned n_rows, const unsigned n_cols,
	dbj_matrix_data_type a[static n_rows][n_cols], dbj_matrix_data_type m[static n_cols][n_rows])
{
	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j)
			m[j][i] = a[i][j];
	return m;
}

DBJ_API dbj_matrix_data_type simple_sdot_sse
(int n, const dbj_matrix_data_type x[static n], const dbj_matrix_data_type y[static n])
{
	int i, n8 = n >> 3 << 3;
	__m128 vs1, vs2;
	dbj_matrix_data_type s, t[4];
	vs1 = _mm_setzero_ps();
	vs2 = _mm_setzero_ps();
	for (i = 0; i < n8; i += 8)
	{
		__m128 vx1, vx2, vy1, vy2;
		vx1 = _mm_loadu_ps(&x[i]);
		vx2 = _mm_loadu_ps(&x[i + 4]);
		vy1 = _mm_loadu_ps(&y[i]);
		vy2 = _mm_loadu_ps(&y[i + 4]);
		vs1 = _mm_add_ps(vs1, _mm_mul_ps(vx1, vy1));
		vs2 = _mm_add_ps(vs2, _mm_mul_ps(vx2, vy2));
	}
	for (s = 0.0f; i < n; ++i)
		s += x[i] * y[i];
	_mm_storeu_ps(t, vs1);
	s += t[0] + t[1] + t[2] + t[3];
	_mm_storeu_ps(t, vs2);
	s += t[0] + t[1] + t[2] + t[3];
	return s;
}

DBJ_API void* SSE_matmul_simpler(
	const unsigned n_a_rows,
	const unsigned n_a_cols,
	const unsigned n_b_cols,
	dbj_matrix_data_type a[static n_a_rows][n_a_cols],
	dbj_matrix_data_type b[static n_a_rows][n_b_cols],
	dbj_matrix_data_type m[static n_a_rows][n_b_cols])
{
	const unsigned n_b_rows = n_a_cols;

	dbj_matrix_data_type* Temp; ALLOC_WITH_POLICY(Temp, n_b_cols * n_b_rows * sizeof(dbj_matrix_data_type));
	// Temp rows and cols are inverted !
	typedef dbj_matrix_data_type(*matrix)[n_b_rows];
	matrix bT = (matrix)Temp;
	(void)simple_mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = simple_sdot_sse(n_a_cols, a[i], bT[j]);
	free(Temp);
	return m;
}

DBJ_API void* SSE_matmul_better(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	dbj_matrix_data_type a[static n_a_rows][n_a_cols],
	dbj_matrix_data_type b[static n_a_rows][n_b_cols],
	dbj_matrix_data_type m[static n_a_rows][n_b_cols])
{
	const unsigned x = 16, n_b_rows = n_a_cols;

	dbj_matrix_data_type* Temp; ALLOC_WITH_POLICY(Temp, n_b_cols * n_b_rows * sizeof(dbj_matrix_data_type));
	// Temp rows and cols are inverted !
	typedef dbj_matrix_data_type(*matrix)[n_b_rows];
	matrix bT = (matrix)Temp;
	(void)simple_mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; i += x)
	{
		for (unsigned j = 0; j < n_b_cols; j += x)
		{
			unsigned je = n_b_cols < j + x ? n_b_cols : j + x;
			unsigned ie = n_a_rows < i + x ? n_a_rows : i + x;
			for (unsigned ii = i; ii < ie; ++ii)
				for (unsigned jj = j; jj < je; ++jj)
					m[ii][jj] += simple_sdot_sse(n_a_cols, a[ii], bT[jj]);
		}
	}
	free(Temp);
	return m;
}

#endif // __SSE__
/////////////////////////////////////////////////////////////////////////

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
dbj_matrix_data_type* array1d_first(const dbj_matrix_data_type* a, const unsigned rows_a, const unsigned cols_a,
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
dbj_matrix_data_type* array1d_second(const dbj_matrix_data_type* a, const unsigned rows_a, const unsigned cols_a,
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
			a[i * cols_a + j] = (dbj_matrix_data_type)(i * cols_a + j);
		}
	}
}

DBJ_API void* array_ptr_args(
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

DBJ_API void* matrix_args(
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

	printf("\nTesting various matrix multiplication algorithms"
		"\nAll matrices are square, side size is: %d"
		"\nData type is: %s\n\n", DBJ_MATRIX_SIDE_DIMENSION, dbj_matrix_data_type_name);

#ifdef WIN32
	// VT100 ESC
	system(" ");
#endif
}

__attribute__((destructor)) DBJ_API void app_end(void)
{
	free(app_data.a);
	free(app_data.b);
	free(app_data.r);
}
/////////////////////////////////////////////////////////////////////////

UBENCH(matmul, array1d_first) {
	array1d_first(app_data.a, app_data.rows_a, app_data.cols_a, app_data.b, app_data.rows_b, app_data.cols_b, app_data.r);
}

UBENCH(matmul, array1d_second) {
	array1d_second(app_data.a, app_data.rows_a, app_data.cols_a, app_data.b, app_data.rows_b, app_data.cols_b, app_data.r);
}

UBENCH(matmul, array_ptr_args) {
	array_ptr_args(
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

UBENCH(matmul, matrix_args) {
	matrix_args(
		app_data.rows_a, app_data.cols_a, app_data.rows_b, app_data.cols_b,
		app_data.rows_r, app_data.cols_r, (void*)app_data.a,
		(void*)app_data.b, (void*)app_data.r);
}

#if __SSE__

UBENCH(matmul, SSE_simpler) {
	SSE_matmul_simpler(
		app_data.rows_a, app_data.cols_a, app_data.cols_b, (void*)app_data.a,
		(void*)app_data.b, (void*)app_data.r);
}

UBENCH(matmul, SSE_better) {
	SSE_matmul_better(
		app_data.rows_a, app_data.cols_a, app_data.cols_b,
		(void*)app_data.a, (void*)app_data.b, (void*)app_data.r);
}

#endif // __SSE__

/////////////////////////////////////////////////////////////////////////

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
