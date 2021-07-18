/*
 (c) 2021 by dbj at dbj dot org -- https://dbj.org/license_dbj/

https://godbolt.org/z/4zWs9MhP7
https://godbolt.org/z/1KTE3PnEP

 This is benchmarking of a collection of matrix multiplication algorithms.
 Algorithms are kept as simple as possible. No structs are passed as arguments.
 No "clever" "generic" matrix macros are used

 Different compilers multiplied with different platforms multiplied with selection
 of data types  yield a complex picture of benchmarking results.

 Although here is strong hint for you: The simplest algorithm is the fastest.
 Keep in mind compiler has the easiest job optimizing the simplest code.

For smaller matrices it is almost irrelevant which algorithm is used
Were what is "smaller" depends on the runtime.

Generally up to 1024 x 1024 the best algorithm from here would suffice
above that use BLAS or LAPAC or one of the many variants of them two.
GPU's have matmuls "on board" so that is another option; for that see CUDA

 Use this file to recompile and re measure whenever selecting
 the right matrix multiplication algorithm for a new platform

 Here is one run on a WIN10 PRO machine, 8GB RAM , and one very old i5 CPU model

[==========] Running 5 benchmarks.
[ RUN      ] matmul.the_most_by_the_book_matrix_mult
[       OK ] matmul.the_most_by_the_book_matrix_mult (mean 1.885s, confidence interval +- 0.117107%)
[ RUN      ] matmul.matmul_mx_as_array
[       OK ] matmul.matmul_mx_as_array (mean 1.956s, confidence interval +- 0.795940%)
[ RUN      ] matmul.matmul_mx_as_array_another
[       OK ] matmul.matmul_mx_as_array_another (mean 227.829ms, confidence interval +- 0.472907%)
[ RUN      ] matmul.matmul_sdot8
[       OK ] matmul.matmul_sdot8 (mean 413.246ms, confidence interval +- 0.523624%)
[ RUN      ] matmul.matmul_sdot_1
[       OK ] matmul.matmul_sdot_1 (mean 237.796ms, confidence interval +- 0.628446%)
[==========] 5 benchmarks ran.
[  PASSED  ] 5 benchmarks.

OMP. What OMP? OMP might help but on *very* large data sets. Were using GPU is much more reasonable anyway.
Here is a code you can try on your machine with large `DATA_SIZE`: https://godbolt.org/z/fTqTsc5vs


Matrix multiplication required dimensions relationship

A(3,2) x B(2,3) = R(3,3)
							  +------+------+------+
							  |      |      |      |
							  |      |      |      |
							  |      |      |      |
						 B    +--------------------+  2 rows
							  |      |      |      |
							  |      |      |      |
							  |      |      |      |
							  +------+------+------+

			2 columns              3 columns

		+------+------+       +------+------+------+
		|      |      |       |      |      |      |
3 rows  |      |      |       |      |      |      |
		|      |      |       |      |      |      |
		+-------------+       +--------------------+
		|      |      |       |      |      |      |
A       |      |      |   R   |      |      |      |  3 rows
		|      |      |       |      |      |      |
		+-------------+       +--------------------+
		|      |      |       |      |      |      |
		|      |      |       |      |      |      |
		|      |      |       |      |      |      |
		+------+------+       +------+------+------+
*/

// larger side is  * 2
// ignored for testing, for testing see the data used bellow
#define DBJ_MX_SMALLER_SIDE 256

// A * B = R
// few algortihms are using transposed B
// if A and B are not changed between the calls
// another optimisation is to "pre transpose" it once
// Ditto
#define DBJ_MX_ALREADY_TRANSPOSED 1

// for testing make it 0
#define DBJ_BENCHMARKING 1
// when on godbolt make it 1
// godbolt times out for even small sizes
// is is used just to confirm code compiles with no warnings
#define DBJ_ON_GODBOLT 0

#ifdef _MSC_VER
#pragma region common trash
#endif

#define FOR(C, R) for (unsigned C = 0; C < R; ++C)

/* NDEBUG == RELEASE */
#include <assert.h>

#if (defined(__clang__) || defined(__GNUC__))
#define DBJ_CLANGNUC 1
#else
#define DBJ_CLANGNUC 0
#endif

#if !DBJ_ON_GODBOLT
#include "build_time_stamp.inc" // DBJ_BUILD_TIMESTAMP
#if DBJ_BENCHMARKING
#include "ubench.h/ubench.h"
#else
#include "utest.h/utest.h"
#endif // ! DBJ_BENCHMARKING

#else // on godbolt

#if DBJ_BENCHMARKING
#include "https://raw.githubusercontent.com/sheredom/ubench.h/master/ubench.h"
#else
#include "https://raw.githubusercontent.com/sheredom/utest.h/master/utest.h"
#endif // ! DBJ_BENCHMARKING
#define DBJ_BUILD_TIMESTAMP __DATE__ " " __TIME__

#endif // DBJ_ON_GODBOLT

#define DBJ_VT_RESET "\033[0m"
#define DBJ_VT_GREEN "\033[32m"
#define DBJ_VT_RED "\033[31m"

// currently ctor/dtor feature is not used
#if DBJ_CLANGNUC
#define DBJ_CTOR __attribute__((constructor))
#define DBJ_DTOR __attribute__((destructor))
#else
#define DBJ_CTOR
#define DBJ_DTOR
#endif

#ifdef _MSC_VER
#pragma endregion // common trash

#pragma region common data
#endif
//
// dimensions defintions

#if DBJ_BENCHMARKING

// NOTE: Be carefull with sizes
//       UBENCH repeats execution so matrix size is not the prevailing factor
//       keep them small-ish
#define DBJ_MX_A_ROWS DBJ_MX_SMALLER_SIDE
#define DBJ_MX_A_COLS DBJ_MX_SMALLER_SIDE * 2
#define DBJ_MX_B_ROWS DBJ_MX_A_COLS
#define DBJ_MX_B_COLS DBJ_MX_A_ROWS

#else // testing
/*
		 In case of testing we use this constelation
		 of matrices, to check the correctness of algorithms

	 *     ! 1 2 |      | 5 6 |       | 19 22 |
	 *     |     |  x   |     |  =    |       |
	 *     | 3 4 |      | 7 8 |       | 43 50 |
	 */
#define DBJ_MX_A_ROWS 2
#define DBJ_MX_A_COLS 2
#define DBJ_MX_B_ROWS DBJ_MX_A_COLS
#define DBJ_MX_B_COLS 2

#endif // testing

#define DBJ_MX_R_ROWS DBJ_MX_A_ROWS
#define DBJ_MX_R_COLS DBJ_MX_B_COLS

static_assert(DBJ_MX_A_COLS == DBJ_MX_B_ROWS, "DBJ_MX_A_COLS != DBJ_MX_B_ROWS");
static_assert(DBJ_MX_A_ROWS == DBJ_MX_R_ROWS, "DBJ_MX_A_ROWS != DBJ_MX_R_ROWS");
static_assert(DBJ_MX_B_COLS == DBJ_MX_R_COLS, "DBJ_MX_B_COLS != DBJ_MX_R_COLS");

typedef double dbj_matrix_data_type;
#define dbj_matrix_data_type_name "double"

// ubench functions have no parameters
// thus we use common data aka globals
typedef struct app_data_struct
{
	unsigned rows_a;
	unsigned cols_a;
	unsigned rows_b;
	unsigned cols_b;
	unsigned rows_r;
	unsigned cols_r;
	// transposed B dimension
	unsigned rows_bT;
	unsigned cols_bT;
	// the matrixes
	dbj_matrix_data_type a[DBJ_MX_A_ROWS][DBJ_MX_A_COLS];
	dbj_matrix_data_type b[DBJ_MX_B_ROWS][DBJ_MX_B_COLS];
	// transposed b
	dbj_matrix_data_type bT[DBJ_MX_B_COLS][DBJ_MX_B_ROWS];
	// the result
	dbj_matrix_data_type r[DBJ_MX_R_ROWS][DBJ_MX_R_COLS];

} app_data_type;

static inline void* reset_test_result(app_data_type* app_data_)
{
	return memset((void*)app_data_->r, 0, sizeof(dbj_matrix_data_type[DBJ_MX_R_ROWS * DBJ_MX_R_COLS]));
}

#ifdef _MSC_VER
#pragma endregion // common data
#pragma region matrix functions and various matmuls
#endif

#if DBJ_BENCHMARKING

static void* matrix_arr_init(const unsigned rows_a, const unsigned cols_a, dbj_matrix_data_type a[static rows_a][cols_a])
{
	for (unsigned i = 0; i < rows_a; i++)
	{
		for (unsigned j = 0; j < cols_a; j++)
		{
			a[i][j] = (dbj_matrix_data_type)(i * cols_a + j);
		}
	}
	return a;
}
#endif // DBJ_BENCHMARKING

#define dbj_matrix_size_bytes(rows_, cols_, type_) (rows_ * cols_ * sizeof(type_))

static void dbj_matrix_transpose(
	const unsigned rows_m,
	const unsigned cols_m,
	const dbj_matrix_data_type m[static rows_m][cols_m],
	dbj_matrix_data_type t[static cols_m][rows_m])
{
	for (size_t i = 0; i < rows_m; i++)
	{
		for (size_t j = 0; j < cols_m; j++)
		{
			t[j][i] = m[i][j];
		}
	}
}

inline dbj_matrix_data_type sdot_1(int n, const dbj_matrix_data_type x[static n], const dbj_matrix_data_type y[static n])
{
	dbj_matrix_data_type s = (dbj_matrix_data_type)0;
	for (int i = 0; i < n; ++i)
		s += x[i] * y[i];
	return s;
}

inline dbj_matrix_data_type sdot_8(int n, const dbj_matrix_data_type x[static n], const dbj_matrix_data_type y[static n])
{
	int i, n8 = n >> 3 << 3;
	dbj_matrix_data_type s = (dbj_matrix_data_type)0, t[8] = { (dbj_matrix_data_type)0 };
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
	for (s = (dbj_matrix_data_type)0; i < n; ++i)
		s += x[i] * y[i];
	s += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
	return s;
}

// the most "by the book" C matrix mutliplication function
// the key fact might be this is the matrix mutliplication so
// "severley optimized" by compilers there is no point investing
// in finding faster algorithms, including SSE/AVX usage
// for small matrices of course
static void the_most_by_the_book_matrix_mult(
	size_t a_rows,
	size_t a_cols,
	size_t b_cols,
	dbj_matrix_data_type A[static a_rows][a_cols],
	dbj_matrix_data_type B[static a_cols][b_cols],
	dbj_matrix_data_type C[static a_rows][b_cols])
{
	for (size_t i = 0; i < a_rows; ++i)
	{
		for (size_t j = 0; j < b_cols; ++j)
		{
			C[i][j] = 0.0;
			for (size_t l = 0; l < a_cols; ++l)
			{
				C[i][j] += A[i][l] * B[l][j];
			}
		}
	}
}

/*
 use 1D aray as matrix type
 calculate the 1D index of "matrix" [row][col]
 this is in here because it is curiously and persistently the fastest matmul
 when ad if no prepocessing is allowed
 and when matrixes are small-ish
 */
static dbj_matrix_data_type* matmul_mx_as_array(
	const size_t a_rows, const size_t a_cols, const size_t b_cols,
	dbj_matrix_data_type* a, dbj_matrix_data_type* b, dbj_matrix_data_type* c)
{
	/*
	the matmul dimensional requirements
	A rows    == B columns
	A columns == B rows
	R rows    == A rows
	R columns == B columns
	*/
	for (size_t i = 0; i < a_rows; i++)
	{
		for (size_t k = 0; k < b_cols; k++)
		{
			dbj_matrix_data_type sum = (dbj_matrix_data_type)0.0;
			for (size_t j = 0; j < a_cols /* same as b rows */; j++)
			{
				sum += a[i * a_cols + j] * b[j * a_rows + k];
			}
			c[i * a_rows + k] = sum;
		}
	}
	return c;
}

/* ---------------------------------------------------------------------------- */
static dbj_matrix_data_type* matmul_mx_as_array_another(
	const size_t a_rows,
	const size_t a_cols,
	const size_t b_cols,
	dbj_matrix_data_type* a,
	dbj_matrix_data_type* b,
	dbj_matrix_data_type* c,
	dbj_matrix_data_type* bT)
{
	// orinteering
	// const unsigned b_rows  = a_cols;
	// const unsigned bt_rows = b_cols;
	// const unsigned bt_cols = b_rows ;

	dbj_matrix_data_type* bTR = bT;
#if DBJ_MX_ALREADY_TRANSPOSED == 0
	dbj_matrix_transpose(a_cols, b_cols, (void*)b, (void*)bTR);
#endif

	for (unsigned i = 0; i < a_rows; i++)
	{
		for (unsigned k = 0; k < b_cols; k++)
		{
			dbj_matrix_data_type sum = 0.0;
			for (unsigned j = 0; j < a_cols; j++)
			{
				sum += a[i * a_cols + j] * bTR[k * b_cols + j];
			}
			c[i * b_cols + k] = sum;
		}
	}
	return c;
}

// this is VMT based
static void* matmul_sdot8(
	const unsigned a_rows, const unsigned a_cols, const unsigned b_cols,
	dbj_matrix_data_type a[static a_rows][a_cols],
	dbj_matrix_data_type b[static a_cols][b_cols],
	dbj_matrix_data_type m[static a_rows][b_cols],
	// allocated space for transposed b
	dbj_matrix_data_type bT[static b_cols][a_cols])
{ // orinteering

	dbj_matrix_data_type(*bTR)[a_cols] = bT; // row view on bT

#if DBJ_MX_ALREADY_TRANSPOSED == 0
	dbj_matrix_transpose(a_cols, b_cols, (void*)b, (void*)bTR);
#endif

	for (unsigned i = 0; i < a_rows; ++i)
		for (unsigned j = 0; j < b_cols; ++j)
			m[i][j] = sdot_8(a_cols, a[i], bTR[j]);

	return m;
}

static void* matmul_sdot_1(
	const unsigned a_rows, const unsigned a_cols, const unsigned b_cols,
	dbj_matrix_data_type a[static a_rows][a_cols],
	dbj_matrix_data_type b[static a_cols][b_cols],
	dbj_matrix_data_type m[static a_rows][b_cols],
	// allocated space for transposed b
	dbj_matrix_data_type bT[static b_cols][a_cols])
{ // orienteering
	// const unsigned b_rows  = a_cols;
	// const unsigned bt_rows = b_cols;
	// const unsigned bt_cols = b_rows ;

	// pointer to bT Row
	dbj_matrix_data_type(*bTR)[a_cols] = bT;

#if DBJ_MX_ALREADY_TRANSPOSED == 0
	dbj_matrix_transpose(a_cols, b_cols, (void*)b, (void*)bTR);
#endif

	for (unsigned i = 0; i < a_rows; ++i)
		for (unsigned j = 0; j < b_cols; ++j)
			m[i][j] = sdot_1(a_cols, a[i], bTR[j]);
	return m;
}

/*
https://raw.githubusercontent.com/ShrohanMohapatra/matrix_multiply_quadratic/master/matrix_multiply_test.py

(c) 20201 by dbj@dbj.org -- transformation to C
 */

static void mohapatra(
	const unsigned Arows, /* == Bcols */
	const unsigned Acols, /* == Brows */
	double A[Arows][Acols], double B[Acols][Arows], double E[Arows][Arows])
{
	double maxi = 0;
	FOR(i, Arows)
		FOR(j, Acols)
	{
		if (maxi < A[i][j])
			maxi = A[i][j];
		if (maxi < B[i][j])
			maxi = B[i][j];
	}
	int M = (int)(log10(maxi)) + 1;
	int P = (int)(log10((pow(10, (2 * M)) - 1) * Arows)) + 1;

	double C[Arows]; /* = { 0 }*/
	;
	double D[Acols] /*= { 0 }*/;

	FOR(i, Arows)
	{
		double sum_1 = 0;
		FOR(j, Acols)
			sum_1 = sum_1 * (pow(10, P)) + A[i][j];
		C[i] = sum_1;
	}
	FOR(j, Acols)
	{
		double sum_1 = 0;
		FOR(i, Arows)
		{
			sum_1 = sum_1 * (pow(10, P)) + B[Arows - 1 - i][j];
		}
		D[j] = sum_1;
	}

	FOR(i, Arows)
	{
		FOR(j, Arows)
		{
			E[i][j] = (int)(C[i] * D[j] / (pow(10, (P * (Arows - 1))))) % (int)(pow(10, P));
		}
	}
} // mohapatra

#ifdef _MSC_VER
#pragma endregion // matrix functions and various matmuls
#pragma region common for testing or benchmarking
#endif

static app_data_type* app_data = 0;

static void app_start(void)
{
	app_data = calloc(1, sizeof(app_data_type));
	assert(app_data);

	app_data->rows_a = DBJ_MX_A_ROWS;
	app_data->cols_a = DBJ_MX_A_COLS;
	app_data->rows_b = DBJ_MX_B_ROWS;
	app_data->cols_b = DBJ_MX_B_COLS;
	// transposed b
	app_data->rows_bT = DBJ_MX_B_COLS;
	app_data->cols_bT = DBJ_MX_B_ROWS;
	/* the result */
	app_data->rows_r = DBJ_MX_A_ROWS;
	app_data->cols_r = DBJ_MX_B_COLS;

#if !DBJ_BENCHMARKING
	// testing
	/*
	 *     ! 1 2 |      | 5 6 |       | 19 22 |
	 *     |     |  x   |     |  =    |       |
	 *     | 3 4 |      | 7 8 |       | 43 50 |
	 */
	assert(app_data->rows_a * app_data->cols_a == 4);
	assert(app_data->rows_b * app_data->cols_b == 4);
	assert(app_data->rows_r * app_data->cols_r == 4);

	app_data->a[0][0] = 1;
	app_data->a[0][1] = 2;
	app_data->a[1][0] = 3;
	app_data->a[1][1] = 4;

	app_data->b[0][0] = 5;
	app_data->b[0][1] = 6;
	app_data->b[1][0] = 7;
	app_data->b[1][1] = 8;

#endif // !DBJ_BENCHMARKING

#if DBJ_BENCHMARKING

#define DBJ_APP_KIND "BENCHMARKING"

	matrix_arr_init(app_data->rows_a, app_data->cols_a, app_data->a);
	matrix_arr_init(app_data->rows_b, app_data->cols_b, app_data->b);
#else
#define DBJ_APP_KIND "TESTING"
#endif // ! DBJ_BENCHMARKING

#if DBJ_MX_ALREADY_TRANSPOSED == 1
	dbj_matrix_transpose(app_data->rows_b, app_data->cols_b, (void*)app_data->b, (void*)app_data->bT);
#endif // DBJ_MX_ALREADY_TRANSPOSED

	const float size_a = dbj_matrix_size_bytes(app_data->rows_a, app_data->cols_a, dbj_matrix_data_type) / 1024.0f;
	const float size_b = dbj_matrix_size_bytes(app_data->rows_b, app_data->cols_b, dbj_matrix_data_type) / 1024.0f;
	const float size_bT = dbj_matrix_size_bytes(app_data->rows_bT, app_data->cols_bT, dbj_matrix_data_type) / 1024.0f;
	const float size_r = dbj_matrix_size_bytes(app_data->rows_r, app_data->cols_r, dbj_matrix_data_type) / 1024.0f;

	fprintf(stderr, "\n\n" DBJ_VT_RED " " DBJ_APP_KIND " " DBJ_VT_RESET " various matrix multiplication algorithms"
		"\n(c) 2021 by dbj dot org, https://dbj.org/license_dbj \nTimestamp: %s"
		"\n\nMatrices are\n"
		"\nA :%4d * %4d * sizeof(%s) == %4.2f KB"
		"\nB :%4d * %4d * sizeof(%s) == %4.2f KB"
		"\nbT:%4d * %4d * sizeof(%s) == %4.2f KB"
		"\nR :%4d * %4d * sizeof(%s) == %4.2f KB\n\n" DBJ_VT_RESET,
		DBJ_BUILD_TIMESTAMP,
		app_data->rows_a, app_data->cols_a, dbj_matrix_data_type_name, size_a,
		app_data->rows_b, app_data->cols_b, dbj_matrix_data_type_name, size_b,
		app_data->rows_bT, app_data->cols_bT, dbj_matrix_data_type_name, size_bT,
		app_data->rows_r, app_data->cols_r, dbj_matrix_data_type_name, size_r);
#undef DBJ_APP_KIND
}

static void app_end(void)
{
	if (app_data)
	{
		free(app_data);
		app_data = 0;
	}
}
/////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma endregion // common for testing or benchmarking
#endif

#if DBJ_BENCHMARKING

// rezult reset and checking are done in UTEST's, see bellow

UBENCH(matmul, mohapatra)
{
	mohapatra(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS,
		app_data->a, app_data->b, app_data->r);
}

UBENCH(matmul, transpose_sdot_another)
{
	matmul_sdot_1(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data->a, app_data->b, app_data->r, app_data->bT);
}

UBENCH(matmul, transpose_and_sdot8)
{
	matmul_sdot8(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data->a, app_data->b, app_data->r, app_data->bT);
}

UBENCH(matmul, mx_as_array_another)
{
	matmul_mx_as_array_another(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data->a,
		(void*)app_data->b,
		(void*)app_data->r,
		(void*)app_data->bT);
}

UBENCH(matmul, mx_as_array)
{
	matmul_mx_as_array(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data->a,
		(void*)app_data->b,
		(void*)app_data->r);
}

UBENCH(matmul, the_most_by_the_book_matrix_mult)
{
	the_most_by_the_book_matrix_mult(
		DBJ_MX_A_ROWS,
		DBJ_MX_A_COLS,
		DBJ_MX_B_COLS,
		app_data->a,
		app_data->b,
		app_data->r);
}

#else // testing /////////////////////////////////////////////////////
/*
 *     ! 1 2 |      | 5 6 |       | 19 22 |
 *     |     |  x   |     |  =    |       |
 *     | 3 4 |      | 7 8 |       | 43 50 |
 */
#define check_test_input()                                     \
	do                                                         \
	{                                                          \
		EXPECT_EQ(app_data->a[0][0], (dbj_matrix_data_type)1); \
		EXPECT_EQ(app_data->a[0][1], (dbj_matrix_data_type)2); \
		EXPECT_EQ(app_data->a[1][0], (dbj_matrix_data_type)3); \
		EXPECT_EQ(app_data->a[1][1], (dbj_matrix_data_type)4); \
															   \
		EXPECT_EQ(app_data->b[0][0], (dbj_matrix_data_type)5); \
		EXPECT_EQ(app_data->b[0][1], (dbj_matrix_data_type)6); \
		EXPECT_EQ(app_data->b[1][0], (dbj_matrix_data_type)7); \
		EXPECT_EQ(app_data->b[1][1], (dbj_matrix_data_type)8); \
	} while (0)

#define check_test_result()                                     \
	do                                                          \
	{                                                           \
		EXPECT_EQ(app_data->r[0][0], (dbj_matrix_data_type)19); \
		EXPECT_EQ(app_data->r[0][1], (dbj_matrix_data_type)22); \
		EXPECT_EQ(app_data->r[1][0], (dbj_matrix_data_type)43); \
		EXPECT_EQ(app_data->r[1][1], (dbj_matrix_data_type)50); \
	} while (0)

UTEST(matmul, mohapatra)
{
	reset_test_result(app_data);
	mohapatra(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS,
		app_data->a, app_data->b, app_data->r);
	check_test_result();
}

UTEST(matmul, transpose_sdot_another)
{
	reset_test_result(app_data);
	matmul_sdot_1(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data->a, app_data->b, app_data->r, app_data->bT);
	check_test_result();
}

UTEST(matmul, transpose_and_sdot8)
{
	reset_test_result(app_data);
	matmul_sdot8(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data->a, app_data->b, app_data->r, app_data->bT);
	check_test_result();
}

UTEST(matmul, mx_as_array_another)
{
	reset_test_result(app_data);
	matmul_mx_as_array_another(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data->a,
		(void*)app_data->b,
		(void*)app_data->r,
		(void*)app_data->bT);
	check_test_result();
}

UTEST(matmul, mx_as_array)
{
	reset_test_result(app_data);
	matmul_mx_as_array(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data->a,
		(void*)app_data->b,
		(void*)app_data->r);
	check_test_result();
}

UTEST(matmul, the_most_by_the_book_matrix_mult)
{
	reset_test_result(app_data);
	the_most_by_the_book_matrix_mult(
		DBJ_MX_A_ROWS,
		DBJ_MX_B_ROWS,
		DBJ_MX_A_COLS,
		app_data->a,
		app_data->b,
		app_data->r);
	check_test_result();
}

#endif // testing

#ifdef _MSC_VER
#pragma region common main
#endif

#if DBJ_BENCHMARKING
UBENCH_STATE();
#else  // ! DBJ_BENCHMARKING
UTEST_STATE();
#endif // ! DBJ_BENCHMARKING

int main(int argc, const char* const argv[])
{
#if defined(_WIN32)
	// VT100 ESC codes kick-start
	system(" ");
#endif
	app_start();

#if DBJ_BENCHMARKING
	ubench_main(argc, argv);
#else  // ! DBJ_BENCHMARKING
	utest_main(argc, argv);
#endif // ! DBJ_BENCHMARKING

	app_end();

#if defined(_WIN32)
	printf(" " DBJ_VT_RESET " ");
#endif
	return 42;
}
#ifdef _MSC_VER
#pragma endregion // common main
#endif

// EOF
