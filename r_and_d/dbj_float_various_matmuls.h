#ifndef DBJ_FLOAT_VARIOUS_MATMULS_INC
#define DBJ_FLOAT_VARIOUS_MATMULS_INC

#include "dbj_matrix.h"

#include <stdlib.h>
#include <intrin.h> // #define __SSE__ 1
#ifdef HAVE_CBLAS
#include <cblas.h>
#endif

#ifndef DBJ_API 
#define DBJ_API static
#endif // DBJ_API 

//////////////////////////////////////////////////////////////////////
// all of the bellow uses:
typedef DBJ_MATRIX_STRUCT(float) float_matrix_struct;

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

DBJ_API float_matrix_struct* new_float_matrix(const unsigned /*n_rows*/, const unsigned  /*n_cols*/);
DBJ_API float_matrix_struct* make_random_float_matrix(
	const unsigned /*n_rows*/, const unsigned  /*n_cols*/, float (*)(void)
);
DBJ_API void* mat_transpose(
	const unsigned /*n_rows*/, const unsigned /*n_cols*/,
	float /*a*/[/*static n_rows*/][*/*n_cols*/], float /*m*/[/*static n_cols*/][*/*n_rows*/]
);


DBJ_API void* mat_mul_null(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
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

#if __SSE__

DBJ_API void* mat_mul_2(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);

DBJ_API void* mat_mul_7(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);

#endif // __SSE__

#ifdef HAVE_CBLAS
DBJ_API void* mat_mul_5(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);

DBJ_API void* mat_mul_6(const unsigned /*n_a_rows*/, const unsigned /*n_a_cols*/, const unsigned /*n_b_cols*/,
	float a[ /*n_a_rows*/][* /*n_a_cols*/], float b[ /*n_a_rows*/][* /*n_b_cols*/], float m[ /*n_a_rows*/][* /*n_b_cols*/]
);
#endif // HAVE_CBLAS


const description_function_pair algo_table[] =
{
DF_PAIR("0: naive - no optimization" , mat_mul_0) ,
DF_PAIR("1: transposing the second matrix" , mat_mul_1) ,
#if __SSE__
DF_PAIR("2: explicitly vectorized sdot() with SSE" , mat_mul_2) ,
#else
DF_PAIR("2: explicitly vectorized sdot() with SSE not implemented" , mat_mul_null) ,
#endif // ! __SSE__
DF_PAIR("3: implicitly vectorized sdot()" , mat_mul_3) ,
DF_PAIR("4: no vectorization hints", mat_mul_4) ,
#ifdef HAVE_CBLAS
DF_PAIR("5: with sdot() from an external CBLAS library" , mat_mul_5) ,
DF_PAIR("6: with sgemm() from an external CBLAS library", mat_mul_6) ,
#else
DF_PAIR("5: with sdot() from CBLAS library not implemented " , mat_mul_null) ,
DF_PAIR("6: with sgemm() from CBLAS library not implemented ", mat_mul_null) ,
#endif // ! HAVE_CBLAS
#if __SSE__
DF_PAIR("7: explicitly SSE sdot() plus loop tiling", mat_mul_7)
#else
DF_PAIR("7: explicitly SSE sdot() plus loop tiling not implemented" , mat_mul_null) ,
#endif // ! __SSE__
};

//////////////////////////////////////////////////////////////////////
#ifdef DBJ_FLOAT_VARIOUS_MATMULS_IMPLEMENTATION
//////////////////////////////////////////////////////////////////////


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
DBJ_API float_matrix_struct* make_random_float_matrix(
	const unsigned n_rows, const unsigned  n_cols, float (*float_rand)(void)
)
{
	float_matrix_struct* fmt_pointer = new_float_matrix(n_rows, n_cols);

	DBJ_MATRIX_ALIAS(matrix, float, n_cols);
	DBJ_MATRIX_CAST(mx, matrix, fmt_pointer);

	for (unsigned i = 0; i < n_rows; ++i)
		for (unsigned j = 0; j < n_cols; ++j)
			mx[i][j] = float_rand(); //
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

#ifdef __SSE__

DBJ_API void* mat_mul_2(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
)
{
	const unsigned n_b_rows = n_a_cols;

	float_matrix_struct* Temp = new_float_matrix(n_b_cols, n_b_rows);
	DBJ_MATRIX_ALIAS(matrix, float, n_b_rows);
	DBJ_MATRIX_CAST(bT, matrix, Temp);
	(void)mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_sse(n_a_cols, a[i], bT[j]);
	free(Temp);
	return m;
}

DBJ_API void* mat_mul_7(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
)
{
	const unsigned x = 16, n_b_rows = n_a_cols;

	float_matrix_struct* Temp = new_float_matrix(n_b_cols, n_b_rows);
	DBJ_MATRIX_ALIAS(matrix, float, n_b_rows);
	DBJ_MATRIX_CAST(bT, matrix, Temp);
	(void)mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; i += x) {
		for (unsigned j = 0; j < n_b_cols; j += x) {
			int je = n_b_cols < j + x ? n_b_cols : j + x;
			int ie = n_a_rows < i + x ? n_a_rows : i + x;
			for (unsigned ii = i; ii < ie; ++ii)
				for (unsigned jj = j; jj < je; ++jj)
					m[ii][jj] += sdot_sse(n_a_cols, a[ii], bT[jj]);
		}
	}
	free(Temp);
	return m;
}
#endif // __SSE__

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

DBJ_API void* mat_mul_5(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
)
{
	const unsigned n_b_rows = n_a_cols;

	float_matrix_struct* Temp = new_float_matrix(n_b_cols, n_b_rows);
	DBJ_MATRIX_ALIAS(matrix, float, n_b_rows);
	DBJ_MATRIX_CAST(bT, matrix, Temp);
	(void)mat_transpose(n_b_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = cblas_sdot(n_a_cols, a[i], 1, bT[j], 1); // clas_sdot() ???
	free(Temp);
	return m;
}

DBJ_API void* mat_mul_6(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	float a[static n_a_rows][n_a_cols], float b[static n_a_rows][n_b_cols], float m[static n_a_rows][n_b_cols]
)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_a_rows, n_b_cols, n_a_cols, 1.0f, a[0], n_a_rows, b[0], n_b_rows, 0.0f, m[0], n_a_rows);
	return m;
}
#endif // HAVE_CBLAS

//////////////////////////////////////////////////////////////////////
#endif // DBJ_FLOAT_VARIOUS_MATMULS_IMPLEMENTATION
//////////////////////////////////////////////////////////////////////

#undef DF_PAIR

#endif //  DBJ_FLOAT_VARIOUS_MATMULS_INC