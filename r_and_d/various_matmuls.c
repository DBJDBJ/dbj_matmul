#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <intrin.h> // #define __SSE__ 1

#ifdef HAVE_CBLAS
#include <cblas.h>
#endif

#include "pseudo_random.h"
#include "dbj_matrix.h"

#define DBJ_API static

//////////////////////////////////////////////////////////////////////

/*
* 1 sec == 1000 milli sec == 1e+6 micro sec == 1e+9 nano sec
*/

typedef struct { unsigned long v; } dbj_milsec;
typedef struct { unsigned long v; } dbj_sec;

#define dbj_one_sec (dbj_milsec){ .v = 1000 }

#ifdef _WIN32
// from synchapi.h
extern void Sleep(
	unsigned long /*DWORD*/ dwMilliseconds
);
#else
#include <unistd.h>
#endif

DBJ_API void dbj_sleep(dbj_milsec milisec_)
{
#ifdef _WIN32
	Sleep(milisec_.v);
#else
	usleep(milisec_.v * 1000);  /* sleep for 100 milliSeconds */
#endif
}

//////////////////////////////////////////////////////////////////////

typedef void* (*mat_mul_function) (
	const unsigned, const unsigned, const unsigned,
	float a[][*], float b[][*], float m[][*]
	);

typedef struct {
	const char* description;
	mat_mul_function function;
} description_function_pair;

#define DF_PAIR(D_,F_) (description_function_pair){ .description = D_ , .function = F_ }

//////////////////////////////////////////////////////////////////////

DBJ_API void* mat_mul_null(
	const unsigned n_a_rows,
	const unsigned n_a_cols,
	const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols],
	float b[static n_a_rows][n_b_cols],
	float m[static n_a_rows][n_b_cols]
) {
	// used for when algorithms are not implemented yet
	(void)a;	(void)b;	(void)m;	return NULL;
}

DBJ_API void* mat_mul_0(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);

DBJ_API void* mat_mul_1(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);

DBJ_API void* mat_mul_3(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);

DBJ_API void* mat_mul_4(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);

const description_function_pair algo_table[] =
{
DF_PAIR("0: naive - no optimization" , mat_mul_0) ,
DF_PAIR("1: transposing the second matrix" , mat_mul_1) ,
#if __SSE__
DF_PAIR("2: explicitly vectorized sdot() with SSE" , mat_mul_0) ,
#else
DF_PAIR("2: explicitly vectorized sdot() with SSE not implemented" , mat_mul_null) ,
#endif // ! __SSE__
DF_PAIR("3: implicitly vectorized sdot()" , mat_mul_3) ,
DF_PAIR("4: no vectorization hints", mat_mul_4) ,
#ifdef HAVE_CBLAS
DF_PAIR("5: with sdot() from an external CBLAS library" , mat_mul_0) ,
DF_PAIR("6: with sgemm() from an external CBLAS library", mat_mul_0) ,
#else
DF_PAIR("5: with sdot() from CBLAS library not implemented " , mat_mul_null) ,
DF_PAIR("6: with sgemm() from CBLAS library not implemented ", mat_mul_null) ,
#endif // ! HAVE_CBLAS
#if __SSE__
DF_PAIR("7: explicitly SSE sdot() plus loop tiling", mat_mul_0)
#else
DF_PAIR("7: explicitly SSE sdot() plus loop tiling not implemented" , mat_mul_null) ,
#endif // ! __SSE__
};

//////////////////////////////////////////////////////////////////////

typedef DBJ_MATRIX_STRUCT(float) float_matrix_struct;

DBJ_API float_matrix_struct* new_float_matrix(const unsigned n_rows, const unsigned  n_cols)
{
	// float_matrix_struct* retval = NULL; // must be init to NULL

	float_matrix_struct* retval = DBJ_MATRIX_ALLOC(n_rows, n_cols, DBJ_MATRIX_STRUCT_SIZE(float, n_rows, n_cols));

	if (retval) {
		retval->rows = n_rows;
		retval->cols = n_cols;
	}
	else {
		fprintf(stderr, "\n%s(%d) DBJ_MATRIX_NEW failed?\n", __FILE__, __LINE__);
		perror(" ");
		exit(1);
	}
	return retval;
}

/*
allocate new matrix and populate it with random float's
*/
DBJ_API float_matrix_struct* make_random_float_matrix(const unsigned n_rows, const unsigned  n_cols)
{
	float_matrix_struct* fmt_pointer = new_float_matrix(n_rows, n_cols);

	DBJ_MATRIX_ALIAS(matrix, float, n_cols);
	DBJ_MATRIX_CAST(mx, matrix, fmt_pointer);

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
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
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

DBJ_API void* mat_mul_1(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
)
{
	const unsigned n_b_rows = n_a_cols;
	// Temp rows and cols are inverted
	float_matrix_struct* Temp = new_float_matrix(n_b_cols, n_b_rows);
	DBJ_MATRIX_ALIAS(matrix, float, n_b_rows);
	DBJ_MATRIX_CAST(bT, matrix, Temp);
	(void)mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i) {
		const float* ai = a[i];
		float* mi = m[i];
		for (unsigned j = 0; j < n_b_cols; ++j) {
			float t = 0.0f, * bTj = bT[j];
			for (unsigned k = 0; k < n_a_cols; ++k)
				t += ai[k] * bTj[k];
			mi[j] = t;
		}
	}

	free(Temp);

	return m;
}

/////////////////////////////////////////////////////////////////////////////////////
#if 0
/////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////
#endif // 0
/////////////////////////////////////////////////////////////////////////////////////

DBJ_API void* mat_mul_3(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
)
{
	int  n_b_rows = n_a_cols;

	float_matrix_struct* Temp = new_float_matrix(n_b_cols, n_b_rows);
	DBJ_MATRIX_ALIAS(matrix, float, n_b_rows);
	DBJ_MATRIX_CAST(bT, matrix, Temp);
	(void)mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_8(n_a_cols, a[i], bT[j]);
	free(Temp);
	return m;
}

DBJ_API void* mat_mul_4(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
)
{
	int n_b_rows = n_a_cols;

	float_matrix_struct* Temp = new_float_matrix(n_b_cols, n_b_rows);
	DBJ_MATRIX_ALIAS(matrix, float, n_b_rows);
	DBJ_MATRIX_CAST(bT, matrix, Temp);
	(void)mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_1(n_a_cols, a[i], bT[j]);
	free(Temp);
	return m;
}

#ifdef HAVE_CBLAS

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
#endif // HAVE_CBLAS

/**********************************************************************************/

enum {
	default_multiplication_algorithm = 4, max_matrix_side = DBJ_SANITY_MAX_ROW_COUNT
};

/*
These are square matrices, matrix_side_length is the length of the side of that square
matrix data type is pointer to pointer
matrix element type is float
*/
DBJ_API int test_matmul(unsigned matrix_side_length, unsigned algorithm_id, description_function_pair dfp)
{
	int algo = algorithm_id; (void)algo;
	const char* algo_name = dfp.description;
	mat_mul_function mm_fun = dfp.function;

	assert(algo_name);

	DBJ_MATRIX_ALIAS(matrix, float, matrix_side_length);

	float_matrix_struct* mx_struct_A = NULL, * mx_struct_B = NULL, * mx_struct_M = NULL;

	assert(matrix_side_length < max_matrix_side);
	// release mode: quietly adjusting to max side .. not exactly good but ok in this context
	matrix_side_length = matrix_side_length % (size_t)max_matrix_side;

	mx_struct_A = make_random_float_matrix(matrix_side_length, matrix_side_length);
	mx_struct_B = make_random_float_matrix(matrix_side_length, matrix_side_length);

	mx_struct_M = new_float_matrix(matrix_side_length, matrix_side_length);

	DBJ_MATRIX_CAST(mx_A, matrix, mx_struct_A);
	DBJ_MATRIX_CAST(mx_B, matrix, mx_struct_B);
	DBJ_MATRIX_CAST(mx_M, matrix, mx_struct_M);


	clock_t start_time_ = clock();

	(void)mm_fun(
		matrix_side_length, matrix_side_length, matrix_side_length,
		mx_A, mx_B, mx_M
	);

	///////////////////////////////////////////////
#if 0
///////////////////////////////////////////////
#ifdef __SSE__
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
	//////////////////////////////////////////////
#endif // 0
///////////////////////////////////////////////

fprintf(stderr, "\nAlgorithm %-36s: ", algo_name);
fprintf(stderr, "\nCPU time: %2.3g sec", (double)(clock() - start_time_) / CLOCKS_PER_SEC);
fprintf(stderr, "\nCentral cell: %g", mx_M[matrix_side_length / 2][matrix_side_length / 2]);

if (mx_struct_M) free(mx_struct_M);
if (mx_struct_A) free(mx_struct_A);
if (mx_struct_B) free(mx_struct_B);

return EXIT_SUCCESS;
}

static void test_various_matmuls(
	const unsigned mx_size /* = 0xFF */,
	const unsigned outer_loop_ /* reserved */
)
{
	(void)outer_loop_;
	fprintf(stderr, "\n------------------------------------------------");
	fprintf(stderr, "\nMatrix width == height == size: %4d", mx_size);

	static const unsigned algo_table_count = sizeof(algo_table) / sizeof(algo_table[0]);

	for (unsigned k = 0; k < algo_table_count; ++k)
	{
		fprintf(stderr, "\n------------------------------------------------");
		test_matmul(mx_size, k, algo_table[k]);

		dbj_sleep(dbj_one_sec);
#if 0
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
#endif // 0
	}

}

int various_matmuls(const int argc, const char** argv)
{
	(void)argc;
	(void)argv;
	//float_matrix_struct* fmtp = make_random_float_matrix(3, 4);
	//DBJ_MATRIX_FREE(fmtp);
	test_various_matmuls(DBJ_SANITY_MAX_ROW_COUNT / 100, 0 /* reserved */);
	return 42;
}