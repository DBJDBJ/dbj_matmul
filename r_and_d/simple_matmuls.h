#ifndef DBJ_VARIOUS_SIMPLE_MATMULS_INC
#define DBJ_VARIOUS_SIMPLE_MATMULS_INC

// not using any kind of "matrix framework"

#include "../dbj_matmul_common.h" // DBJ_MATRIX_DATA_TYPE

#include <stdlib.h>
#include <intrin.h> // #define __SSE__ 1
#ifdef HAVE_CBLAS
#include <cblas.h>
#endif

#undef DBJ_API
#define DBJ_API

//////////////////////////////////////////////////////////////////////
#ifndef DBJ_SANITY_MAX_ROW_COUNT
#define DBJ_SANITY_MAX_ROW_COUNT 0xFFFF
#endif

#ifndef DBJ_SANITY_MAX_COL_COUNT
#define DBJ_SANITY_MAX_COL_COUNT 0xFFFF
#endif
//////////////////////////////////////////////////////////////////////

typedef void* (*simple_mat_mul_function)(
	const unsigned, const unsigned, const unsigned,
	DBJ_MATRIX_DATA_TYPE[][*],
	DBJ_MATRIX_DATA_TYPE[][*],
	DBJ_MATRIX_DATA_TYPE[][*] /* the result */
	);

typedef struct
{
	const char* description;
	simple_mat_mul_function function;
} simple_description_function_pair;

#define DF_PAIR(D_, F_) \
	(simple_description_function_pair) { .description = D_, .function = F_ }

//////////////////////////////////////////////////////////////////////

DBJ_API DBJ_MATRIX_DATA_TYPE* new_simple_matrix(const unsigned, const unsigned);

DBJ_API DBJ_MATRIX_DATA_TYPE* make_simple_matrix(
	const unsigned, const unsigned, DBJ_MATRIX_DATA_TYPE(*)(void));

DBJ_API void* simple_mat_transpose(
	const unsigned, const unsigned, DBJ_MATRIX_DATA_TYPE[][*], DBJ_MATRIX_DATA_TYPE[][*]);

DBJ_API void* simple_mat_mul_null(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
	// used for when algorithms are not implemented yet
	(void)a;
	(void)b;
	(void)m;
	return NULL;
}

DBJ_API void* simple_mat_mul_0_0(
	const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE*, DBJ_MATRIX_DATA_TYPE*, DBJ_MATRIX_DATA_TYPE*);

DBJ_API void* simple_mat_mul_0(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/]);

DBJ_API void* simple_mat_mul_1(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/]);

DBJ_API void* simple_mat_mul_3(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/]);

DBJ_API void* simple_mat_mul_4(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE[/*n_a_rows*/][* /*n_b_cols*/]);

#if __SSE__

DBJ_API void* simple_mat_mul_2(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE a[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE b[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE m[/*n_a_rows*/][* /*n_b_cols*/]);

DBJ_API void* simple_mat_mul_7(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE a[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE b[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE m[/*n_a_rows*/][* /*n_b_cols*/]);

#endif // __SSE__

#ifdef HAVE_CBLAS
DBJ_API void* simple_mat_mul_5(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE a[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE b[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE m[/*n_a_rows*/][* /*n_b_cols*/]);

DBJ_API void* simple_mat_mul_6(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	DBJ_MATRIX_DATA_TYPE a[/*n_a_rows*/][* /*n_a_cols*/], DBJ_MATRIX_DATA_TYPE b[/*n_a_rows*/][* /*n_b_cols*/], DBJ_MATRIX_DATA_TYPE m[/*n_a_rows*/][* /*n_b_cols*/]);
#endif // HAVE_CBLAS

const simple_description_function_pair dbj_simple_matmuls_algo_table[] =
{
	DF_PAIR("0: MAX optimization", simple_mat_mul_0_0),
	DF_PAIR("1: simple naive - no optimization", simple_mat_mul_0),
	DF_PAIR("2: simple transposing the second matrix", simple_mat_mul_1),
	DF_PAIR("3: simple explicitly vectorized sdot() with SSE", simple_mat_mul_2), /* requires __SSE__ */
	DF_PAIR("4: simple explicitly SSE sdot() plus loop tiling", simple_mat_mul_7) , /* requires __SSE__ */
	DF_PAIR("5: simple implicitly vectorized sdot()", simple_mat_mul_3),
	DF_PAIR("6: simple no vectorization hints", simple_mat_mul_4),
#ifdef HAVE_CBLAS
	DF_PAIR("7: simple with sdot() from an external CBLAS library", simple_mat_mul_5),
	DF_PAIR("8: simple with sgemm() from an external CBLAS library", simple_mat_mul_6),
#endif																	// ! HAVE_CBLAS
};

static const unsigned dbj_simple_matmuls_algo_table_size = (sizeof(dbj_simple_matmuls_algo_table) / sizeof(dbj_simple_matmuls_algo_table[0]));

//////////////////////////////////////////////////////////////////////
#ifdef DBJ_VARIOUS_SIMPLE_MATMULS_IMPLEMENTATION
//////////////////////////////////////////////////////////////////////

DBJ_API DBJ_MATRIX_DATA_TYPE* new_simple_matrix(const unsigned n_rows, const unsigned n_cols)
{
	DBJ_MATRIX_DATA_TYPE* retval = calloc(n_rows * n_cols, sizeof(DBJ_MATRIX_DATA_TYPE));

	if (!retval)
	{
		fprintf(stderr, "\n%s(%d) new_simple_matrix() failed?\n", __FILE__, __LINE__);
		perror(" ");
		exit(1);
	}
	return retval;
}

DBJ_API DBJ_MATRIX_DATA_TYPE* make_simple_matrix(
	const unsigned n_rows, const unsigned n_cols, DBJ_MATRIX_DATA_TYPE(*value_provider)(void))
{
	DBJ_MATRIX_DATA_TYPE* rezult = new_simple_matrix(n_rows, n_cols);

	typedef DBJ_MATRIX_DATA_TYPE(*matrix)[n_cols];
	matrix mx = (matrix)rezult;

	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j)
			mx[i][j] = value_provider(); //
	return rezult;
}

DBJ_API void* simple_mat_transpose(
	const unsigned n_rows, const unsigned n_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_rows][n_cols], DBJ_MATRIX_DATA_TYPE m[static n_cols][n_rows])
{
	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j)
			m[j][i] = a[i][j];
	return m;
}

DBJ_API DBJ_MATRIX_DATA_TYPE simple_sdot_1(int n, const DBJ_MATRIX_DATA_TYPE x[static n], const DBJ_MATRIX_DATA_TYPE y[static n])
{
	DBJ_MATRIX_DATA_TYPE s = 0.0f;
	for (int i = 0; i < n; ++i)
		s += x[i] * y[i];
	return s;
}

DBJ_API DBJ_MATRIX_DATA_TYPE simple_sdot_8(int n, const DBJ_MATRIX_DATA_TYPE x[static n], const DBJ_MATRIX_DATA_TYPE y[static n])
{
	int i, n8 = n >> 3 << 3;
	DBJ_MATRIX_DATA_TYPE s = 0.0f, t[8] = { 0.0f };
	// t[0] = t[1] = t[2] = t[3] = t[4] = t[5] = t[6] = t[7] = 0.0f;
	for (i = 0; i < n8; i += 8)
	{
		t[0] += x[i + 0] * y[i + 0];
		t[1] += x[i + 1] * y[i + 1];
		t[2] += x[i + 2] * y[i + 2];
		t[3] += x[i + 3] * y[i + 3];
		t[4] += x[i + 4] * y[i + 4];
		t[5] += x[i + 5] * y[i + 5];
		t[6] += x[i + 6] * y[i + 6];
		t[7] += x[i + 7] * y[i + 7];
	}
	for (s = 0.0f; i < n; ++i)
		s += x[i] * y[i];
	s += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
	return s;
}

#ifdef __SSE__

DBJ_API DBJ_MATRIX_DATA_TYPE simple_sdot_sse(int n, const DBJ_MATRIX_DATA_TYPE x[static n], const DBJ_MATRIX_DATA_TYPE y[static n])
{
	int i, n8 = n >> 3 << 3;
	__m128 vs1, vs2;
	DBJ_MATRIX_DATA_TYPE s, t[4];
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
#endif // __SSE__

/**************************************************

 Various Matrix multiplication algorithms

 NOTE: check and choose for your RT Env
	   there are differences but are *very* dependant on
	   compiler, OS and hardware

 */

DBJ_API void* simple_mat_mul_0_0(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE* a, DBJ_MATRIX_DATA_TYPE* b, DBJ_MATRIX_DATA_TYPE* m)
{
	// runtime casting takes time
	const unsigned n_b_rows = n_a_cols;
	DBJ_MATRIX_DATA_TYPE(*ax)[n_a_cols] = a;
	DBJ_MATRIX_DATA_TYPE(*bx)[n_b_cols] = b;
	DBJ_MATRIX_DATA_TYPE(*mx)[n_b_rows] = m;

	for (unsigned i = 0; i < n_a_rows; ++i)
	{
		for (unsigned j = 0; j < n_b_cols; ++j)
		{
			DBJ_MATRIX_DATA_TYPE t = 0.0;
			for (unsigned k = 0; k < n_a_cols; ++k)
				t += ax[i][k] * bx[k][j];
			mx[i][j] = t;
		}
	}
	return m;
}

DBJ_API void* simple_mat_mul_0(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
	for (unsigned i = 0; i < n_a_rows; ++i)
	{
		for (unsigned j = 0; j < n_b_cols; ++j)
		{
			DBJ_MATRIX_DATA_TYPE t = 0.0;
			for (unsigned k = 0; k < n_a_cols; ++k)
				t += a[i][k] * b[k][j];
			m[i][j] = t;
		}
	}
	return m;
}

DBJ_API void* simple_mat_mul_1(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
	const unsigned n_b_rows = n_a_cols;

	DBJ_MATRIX_DATA_TYPE* Temp = new_simple_matrix(n_b_cols, n_b_rows);
	// Temp rows and cols are inverted !
	typedef DBJ_MATRIX_DATA_TYPE(*matrix)[n_b_rows];
	matrix bT = (matrix)Temp;
	(void)simple_mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
	{
		const DBJ_MATRIX_DATA_TYPE* ai = a[i];
		DBJ_MATRIX_DATA_TYPE* mi = m[i];
		for (unsigned j = 0; j < n_b_cols; ++j)
		{
			DBJ_MATRIX_DATA_TYPE t = 0.0f, * bTj = bT[j];
			for (unsigned k = 0; k < n_a_cols; ++k)
				t += ai[k] * bTj[k];
			mi[j] = t;
		}
	}

	free(Temp);

	return m;
}

DBJ_API void* simple_mat_mul_2(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
#ifndef __SSE__
#error __SEE__ is required here
#endif

	const unsigned n_b_rows = n_a_cols;

	DBJ_MATRIX_DATA_TYPE* Temp = new_simple_matrix(n_b_cols, n_b_rows);
	// Temp rows and cols are inverted !
	typedef DBJ_MATRIX_DATA_TYPE(*matrix)[n_b_rows];
	matrix bT = (matrix)Temp;
	(void)simple_mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = simple_sdot_sse(n_a_cols, a[i], bT[j]);
	free(Temp);
	return m;
}

DBJ_API void* simple_mat_mul_7(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
#ifndef __SSE__
#error __SEE__ is required here
#endif

	const unsigned x = 16, n_b_rows = n_a_cols;

	DBJ_MATRIX_DATA_TYPE* Temp = new_simple_matrix(n_b_cols, n_b_rows);
	// Temp rows and cols are inverted !
	typedef DBJ_MATRIX_DATA_TYPE(*matrix)[n_b_rows];
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

/////////////////////////////////////////////////////////////////////////////////////

DBJ_API void* simple_mat_mul_3(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
	int n_b_rows = n_a_cols;

	DBJ_MATRIX_DATA_TYPE* Temp = new_simple_matrix(n_b_cols, n_b_rows);
	// Temp rows and cols are inverted !
	typedef DBJ_MATRIX_DATA_TYPE(*matrix)[n_b_rows];
	matrix bT = (matrix)Temp;
	(void)simple_mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = simple_sdot_8(n_a_cols, a[i], bT[j]);
	free(Temp);
	return m;
}

DBJ_API void* simple_mat_mul_4(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
	int n_b_rows = n_a_cols;

	DBJ_MATRIX_DATA_TYPE* Temp = new_simple_matrix(n_b_cols, n_b_rows);
	// Temp rows and cols are inverted !
	typedef DBJ_MATRIX_DATA_TYPE(*matrix)[n_b_rows];
	matrix bT = (matrix)Temp;
	(void)simple_mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = simple_sdot_1(n_a_cols, a[i], bT[j]);
	free(Temp);
	return m;
}

#ifdef HAVE_CBLAS

DBJ_API void* simple_mat_mul_5(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
	const unsigned n_b_rows = n_a_cols;

	DBJ_MATRIX_DATA_TYPE* Temp = new_simple_matrix(n_rows, n_cols);
	// Temp rows and cols are inverted !
	typedef DBJ_MATRIX_DATA_TYPE(*matrix)[n_rows];
	matrix bT = (matrix)Temp;
	(void)simple_mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = cblas_sdot(n_a_cols, a[i], 1, bT[j], 1); // clas_sdot() ???
	free(Temp);
	return m;
}

DBJ_API void* simple_mat_mul_6(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	DBJ_MATRIX_DATA_TYPE a[static n_a_rows][n_a_cols], DBJ_MATRIX_DATA_TYPE b[static n_a_rows][n_b_cols], DBJ_MATRIX_DATA_TYPE m[static n_a_rows][n_b_cols])
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_a_rows, n_b_cols, n_a_cols, 1.0f, a[0], n_a_rows, b[0], n_b_rows, 0.0f, m[0], n_a_rows);
	return m;
}
#endif // HAVE_CBLAS

//////////////////////////////////////////////////////////////////////
#endif // DBJ_VARIOUS_SIMPLE_MATMULS_IMPLEMENTATION
//////////////////////////////////////////////////////////////////////

#undef DF_PAIR

#endif //  DBJ_VARIOUS_SIMPLE_MATMULS_INC