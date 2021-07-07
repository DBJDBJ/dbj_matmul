#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <intrin.h> // #define __SSE__ 1

#include "pseudo_random.h"


#define DBJ_API static

#define SANITY_MAX_ROW_COUNT 0xFFFF
#define SANITY_MAX_COL_COUNT 0xFFFF



/*
 Helper routines for matrix manipulation
*/

#define MATRIC_HANDLE void *

// usage: typedef MATRIX_STRUCT(char) char_matrix ;
#define MATRIX_STRUCT(T_) \
struct {\
	unsigned rows;\
	unsigned cols;\
	T_ data[];\
}

// #define MATRIX_DATA_TYPE(M_) __typeof__(M_->data)

#define MATRIX_STRUCT_SIZE( T_,R_,C_) \
(sizeof( MATRIX_STRUCT(T_) ) + sizeof( T_[R_ * C_] ))

/*
usage:

char_matrix * char_matrix_struct_pointer = NULL ;

MATRIX_NEW( char_matrix_struct_pointer, rows , cols ) ;

assert(char_matrix_struct_pointer) ;

*/

#define MATRIX_NEW(N_,T_,R_,C_) \
do { \
	assert( R_ < SANITY_MAX_ROW_COUNT); \
	assert( C_ < SANITY_MAX_COL_COUNT); \
	 MATRIX_STRUCT(T_) retval = \
	 calloc(1, MATRIX_STRUCT_SIZE(T_,R_,C_) ) ; \
	 if (retval) { \
		 retval->rows = R_ ; \
		 retval->cols = C_ ; \
	 } ; \
	 N_ = retval ; \
} while(0)

#define MATRIX_FREE(M_) do { assert(M_); free(M_); M_ = NULL; } while(0)

#define MATRIX_DATA_POINTER(T_, M_) (T_(*)[M_->rows][M_->cols])M_->data

// matrix type is used with [][]
#define MATRIX_TYPE(T_, C_) T_(*)[C]
#define MATRIX_ALIAS(N_,T_,C_) \
  typedef T_(*N_)[C_]

// MS_ is matrix struct pointer
// usage: MATRIX_TYPE(mx_struct_pointer) mx = MATRIX_CAST(mx_struct_pointer)
// #define MATRIX_CAST(T_, M_) (T_)MATRIX_DATA_POINTER(T_,M_)
#define MATRIX_CAST(N_,A_, M_) A_ N_ = (A_)(M_->data)

//////////////////////////////////////////////////////////////////////

typedef MATRIX_STRUCT(float) float_matrix_struct;

//////////////////////////////////////////////////////////////////////

DBJ_API float_matrix_struct* new_float_matrix(const unsigned n_rows, const unsigned  n_cols)
{
	assert(n_rows < SANITY_MAX_ROW_COUNT);
	assert(n_cols < SANITY_MAX_COL_COUNT);

	float_matrix_struct* retval = calloc(1, sizeof(float_matrix_struct) + sizeof(float[n_rows * n_cols]));
	assert(retval);
	if (retval) {
		retval->rows = n_rows;
		retval->cols = n_cols;
	};
	return retval;
}

/*
allocate new matrix and populate it with random float's
*/
DBJ_API float_matrix_struct* mat_gen_random(const unsigned n_rows, const unsigned  n_cols)
{
	float_matrix_struct* fmt_pointer = new_float_matrix(n_rows, n_cols);

	//typedef typeof(fmt_pointer->data[0])(*mx_type)[fmt_pointer->cols];

	//mx_type mx = (mx_type)(fmt_pointer->data);

	MATRIX_ALIAS(matrix, float, n_cols);
	MATRIX_CAST(mx, matrix, fmt_pointer);

	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j)
			mx[i][j] = pseudo_rand_float(); //
	return fmt_pointer;
}

DBJ_API void* mat_transpose(
	const unsigned n_rows, const unsigned n_cols,
	float a[static n_rows][n_cols], float m[static n_cols][n_rows]
)
{
	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j)
			m[j][i] = a[i][j];
	return m;
}

DBJ_API float sdot_1(int n, const float x[static n], const float y[static n])
{
	float s = 0.0f;
	for (int i = 0; i < n; ++i) s += x[i] * y[i];
	return s;
}

DBJ_API float sdot_8(int n, const float x[static n], const float y[static n])
{
	int i, n8 = n >> 3 << 3;
	float s = 0.0f, t[8] = { 0.0f };
	// t[0] = t[1] = t[2] = t[3] = t[4] = t[5] = t[6] = t[7] = 0.0f;
	for (i = 0; i < n8; i += 8) {
		t[0] += x[i + 0] * y[i + 0];
		t[1] += x[i + 1] * y[i + 1];
		t[2] += x[i + 2] * y[i + 2];
		t[3] += x[i + 3] * y[i + 3];
		t[4] += x[i + 4] * y[i + 4];
		t[5] += x[i + 5] * y[i + 5];
		t[6] += x[i + 6] * y[i + 6];
		t[7] += x[i + 7] * y[i + 7];
	}
	for (s = 0.0f; i < n; ++i) s += x[i] * y[i];
	s += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
	return s;
}

#ifdef __SSE__


DBJ_API float sdot_sse(int n, const float x[static n], const float y[static n])
{
	int i, n8 = n >> 3 << 3;
	__m128 vs1, vs2;
	float s, t[4];
	vs1 = _mm_setzero_ps();
	vs2 = _mm_setzero_ps();
	for (i = 0; i < n8; i += 8) {
		__m128 vx1, vx2, vy1, vy2;
		vx1 = _mm_loadu_ps(&x[i]);
		vx2 = _mm_loadu_ps(&x[i + 4]);
		vy1 = _mm_loadu_ps(&y[i]);
		vy2 = _mm_loadu_ps(&y[i + 4]);
		vs1 = _mm_add_ps(vs1, _mm_mul_ps(vx1, vy1));
		vs2 = _mm_add_ps(vs2, _mm_mul_ps(vx2, vy2));
	}
	for (s = 0.0f; i < n; ++i) s += x[i] * y[i];
	_mm_storeu_ps(t, vs1);
	s += t[0] + t[1] + t[2] + t[3];
	_mm_storeu_ps(t, vs2);
	s += t[0] + t[1] + t[2] + t[3];
	return s;
}
#endif // __SSE__

/**************************************************

 Various Matrix multiplication algorithms

 NOTE: check and choose for your RT Env
	   there are differences but are *very* dependant on
	   compiler, OS and hardware

 */

DBJ_API void* mat_mul_0(
	const unsigned n_a_rows,
	const unsigned n_a_cols,
	const unsigned n_b_cols,
	float a[n_a_rows][n_a_cols],
	float b[n_a_rows][n_b_cols],
	float m[n_a_rows][n_b_cols]
);

DBJ_API void* mat_mul_0(
	const unsigned n_a_rows,
	const unsigned n_a_cols,
	const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols],
	float b[static n_a_rows][n_b_cols],
	float m[static n_a_rows][n_b_cols]
)
{
	for (unsigned i = 0; i < n_a_rows; ++i) {
		for (unsigned j = 0; j < n_b_cols; ++j) {
			float t = 0.0;
			for (unsigned k = 0; k < n_a_cols; ++k)
				t += a[i][k] * b[k][j];
			m[i][j] = t;
		}
	}
	return m;
}

#if 0

DBJ_API float** mat_mul1(int n_a_rows, int n_a_cols, float* const* a, int n_b_cols, float* const* b)
{
	int i, j, k, n_b_rows = n_a_cols;
	float** m, ** bT;
	m = mat_init(n_a_rows, n_b_cols);
	bT = mat_transpose(n_b_rows, n_b_cols, b);
	for (i = 0; i < n_a_rows; ++i) {
		const float* ai = a[i];
		float* mi = m[i];
		for (j = 0; j < n_b_cols; ++j) {
			float t = 0.0f, * bTj = bT[j];
			for (k = 0; k < n_a_cols; ++k)
				t += ai[k] * bTj[k];
			mi[j] = t;
		}
	}
	mat_destroy(bT);
	return m;
}

#ifdef __SSE__
DBJ_API float** mat_mul2(int n_a_rows, int n_a_cols, float* const* a, int n_b_cols, float* const* b)
{
	int i, j, n_b_rows = n_a_cols;
	float** m, ** bT;
	m = mat_init(n_a_rows, n_b_cols);
	bT = mat_transpose(n_b_rows, n_b_cols, b);
	for (i = 0; i < n_a_rows; ++i)
		for (j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_sse(n_a_cols, a[i], bT[j]);
	mat_destroy(bT);
	return m;
}
DBJ_API float** mat_mul7(int n_a_rows, int n_a_cols, float* const* a, int n_b_cols, float* const* b)
{
	int i, j, ii, jj, x = 16, n_b_rows = n_a_cols;
	float** m, ** bT;
	m = mat_init(n_a_rows, n_b_cols);
	bT = mat_transpose(n_b_rows, n_b_cols, b);
	for (i = 0; i < n_a_rows; i += x) {
		for (j = 0; j < n_b_cols; j += x) {
			int je = n_b_cols < j + x ? n_b_cols : j + x;
			int ie = n_a_rows < i + x ? n_a_rows : i + x;
			for (ii = i; ii < ie; ++ii)
				for (jj = j; jj < je; ++jj)
					m[ii][jj] += sdot_sse(n_a_cols, a[ii], bT[jj]);
		}
	}
	mat_destroy(bT);
	return m;
}
#endif // __SSE__

DBJ_API float** mat_mul3(int n_a_rows, int n_a_cols, float* const* a, int n_b_cols, float* const* b)
{
	int i, j, n_b_rows = n_a_cols;
	float** m, ** bT;
	m = mat_init(n_a_rows, n_b_cols);
	bT = mat_transpose(n_b_rows, n_b_cols, b);
	for (i = 0; i < n_a_rows; ++i)
		for (j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_8(n_a_cols, a[i], bT[j]);
	mat_destroy(bT);
	return m;
}

DBJ_API float** mat_mul4(int n_a_rows, int n_a_cols, float* const* a, int n_b_cols, float* const* b)
{
	int i, j, n_b_rows = n_a_cols;
	float** m, ** bT;
	m = mat_init(n_a_rows, n_b_cols);
	bT = mat_transpose(n_b_rows, n_b_cols, b);
	for (i = 0; i < n_a_rows; ++i)
		for (j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_1(n_a_cols, a[i], bT[j]);
	mat_destroy(bT);
	return m;
}

#ifdef HAVE_CBLAS
#include <cblas.h>

DBJ_API float** mat_mul5(int n_a_rows, int n_a_cols, float* const* a, int n_b_cols, float* const* b)
{
	int i, j, n_b_rows = n_a_cols;
	float** m, ** bT;
	m = mat_init(n_a_rows, n_b_cols);
	bT = mat_transpose(n_b_rows, n_b_cols, b);
	for (i = 0; i < n_a_rows; ++i)
		for (j = 0; j < n_b_cols; ++j)
			m[i][j] = cblas_sdot(n_a_cols, a[i], 1, bT[j], 1);
	mat_destroy(bT);
	return m;
}

DBJ_API** mat_mul6(int n_a_rows, int n_a_cols, float* const* a, int n_b_cols, float* const* b)
{
	float** m, n_b_rows = n_a_cols;
	m = mat_init(n_a_rows, n_b_cols);
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_a_rows, n_b_cols, n_a_cols, 1.0f, a[0], n_a_rows, b[0], n_b_rows, 0.0f, m[0], n_a_rows);
	return m;
}
#endif

/**********************************************************************************/
#endif // 0


enum { default_multiplication_algorithm = 4, max_matrix_side = 0xFF + 0xFF };

/*
These are square matrices, matrix_side_length is the length of the side of that square
matrix data type is pointer to pointer
matrix element type is float
*/
DBJ_API int test_matmul(unsigned matrix_side_length, unsigned algorithm_id, const char* algo_name)
{
	int algo = algorithm_id;

	assert(algo_name);

	MATRIX_ALIAS(matrix, float, matrix_side_length);

	float_matrix_struct* mx_struct_A = NULL, * mx_struct_B = NULL, * mx_struct_M = NULL;

	assert(matrix_side_length < max_matrix_side);
	// release mode: quietly adjusting to max side .. not exactly good but ok in this context
	matrix_side_length = matrix_side_length % (size_t)max_matrix_side;

	mx_struct_A = mat_gen_random(matrix_side_length, matrix_side_length);
	mx_struct_B = mat_gen_random(matrix_side_length, matrix_side_length);

	mx_struct_M = new_float_matrix(matrix_side_length, matrix_side_length);

	MATRIX_CAST(mx_A, matrix, mx_struct_A);
	MATRIX_CAST(mx_B, matrix, mx_struct_B);
	MATRIX_CAST(mx_M, matrix, mx_struct_M);


	clock_t start_time_ = clock();

	if (algo == 0) {
		(void)mat_mul_0(
			matrix_side_length, matrix_side_length, matrix_side_length,
			mx_A, mx_B, mx_M
		);

		// mx_struct_M = mat_mul_0(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
	}
#if 0
	else if (algo == 1) {
		mx_struct_M = mat_mul1(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
#ifdef __SSE__
	}
	else if (algo == 2) {
		mx_struct_M = mat_mul2(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
	}
	else if (algo == 7) {
		mx_struct_M = mat_mul7(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
#endif
	}
	else if (algo == 3) {
		mx_struct_M = mat_mul3(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
	}
	else if (algo == 4) {
		mx_struct_M = mat_mul4(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
#ifdef HAVE_CBLAS
	}
	else if (algo == 5) {
		mx_struct_M = mat_mul5(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
	}
	else if (algo == 6) {
		mx_struct_M = mat_mul6(matrix_side_length, matrix_side_length, mx_struct_A, matrix_side_length, mx_struct_B);
#endif
	}
	else {
		fprintf(stderr, "SKIPPING: unknown algorithm %d\rows_cols", algo);
		goto matmul_exit;
	}
#endif // 0

	fprintf(stderr, "\nAlgorithm %-36s: ", algo_name);
	fprintf(stderr, "CPU time: %2.3g sec", (double)(clock() - start_time_) / CLOCKS_PER_SEC);
	// fprintf(stderr, " Central cell: %g", mx_struct_M[rows_cols / 2][rows_cols / 2]);

	if (mx_struct_M) free(mx_struct_M);
	if (mx_struct_A) free(mx_struct_A);
	if (mx_struct_B) free(mx_struct_B);

	return EXIT_SUCCESS;
}

static const char* algo_name_[] = {
"0: naive - no optimization" ,
"1: transposing the second matrix" ,
"2: explicitly vectorized sdot() with SSE" ,
"3: implicitly vectorized sdot()" ,
"4: no vectorization hints",
"5: with sdot() from an external CBLAS library" ,
"6: with sgemm() from an external CBLAS library",
"7: explicitly SSE sdot() plus loop tiling"
};

static void test_various_matmuls(
	const unsigned mx_size /*= 0xFF*/,
	const unsigned outer_loop_ /*= 0xF*/
)
{
	fprintf(stderr, "\n------------------------------------------------");
	fprintf(stderr, "\nMatrix size: %4d", mx_size);


	for (unsigned k = 0; k < outer_loop_; ++k)
	{
		fprintf(stderr, "\n------------------------------------------------");
		test_matmul(mx_size, 0, algo_name_[0]);
		test_matmul(mx_size, 1, algo_name_[1]);
#ifdef __SSE__
		test_matmul(mx_size, 2, algo_name_[2]);
		test_matmul(mx_size, 7, algo_name_[7]);
#endif
		test_matmul(mx_size, 3, algo_name_[3]);
		test_matmul(mx_size, 4, algo_name_[4]);
#ifdef HAVE_CBLAS
		test_matmul(mx_size, 5, algo_name_[5]);
		test_matmul(mx_size, 6, algo_name_[6]);
#endif 
	}

}

int various_matmuls(const int argc, const char** argv)
{
	(void)argc;
	(void)argv;
	//float_matrix_struct* fmtp = mat_gen_random(3, 4);
	//MATRIX_FREE(fmtp);
	test_various_matmuls(0xF, 0xF);
	return 42;
}