/*

 This is benchmarking of a collection of matrix multiplication algorithms.
 Algorithms are kept as simple as possible. No structs are passed as arguments.
 No "clever" "generic" matrix macros are used

 Different comoilers multiplied with different platforms multiplied selection
 of data types  yield a complex picture of benchmarking results.

 Although here is strong hint for you: The simplest algorithm is the fastest.
 Keep in mind compiler has the easiest job optimizing the simplest code.

 Use this file to recompile and re measure whenever selecting
 the right matrix multiplication algorithm

 https://godbolt.org/z/1KTE3PnEP

 (c) 2021 by dbj at dbj dot org -- https://dbj.org/license_dbj/

 */

#define DBJ_BENCHMARKING 1
#define DBJ_ON_GODBOLT 0


#if ! DBJ_ON_GODBOLT
#include "build_time_stamp.inc" // DBJ_BUILD_TIMESTAMP 
#   if DBJ_BENCHMARKING
#   include "ubench.h/ubench.h"
#else
#   include "utest.h/utest.h"
#endif  // ! DBJ_BENCHMARKING

#else // on godbolt
#   if DBJ_BENCHMARKING
#include "https://raw.githubusercontent.com/sheredom/ubench.h/master/ubench.h"
#else
#   include "https://raw.githubusercontent.com/sheredom/utest.h/master/utest.h"
#endif  // ! DBJ_BENCHMARKING
#define DBJ_BUILD_TIMESTAMP __DATE__ " " __TIME__  
#endif

#pragma region common_trash

#define DBJ_VT_RESET "\033[0m"
#define DBJ_VT_GREEN "\033[32m"
#define DBJ_VT_RED   "\033[31m"

#if (defined(__clang__) || defined(__GNUC__))
#define DBJ_CLANGNUC (1==1)
#else
#define DBJ_CLANGNUC (1==0)
#endif

#if DBJ_CLANGNUC
#define DBJ_CTOR __attribute__((constructor)) 
#define DBJ_DTOR __attribute__((destructor)) 
#else
#define DBJ_CTOR 
#define DBJ_DTOR 
#endif

#define CLANG_IGNORE_PUSH                                                   \
	_Pragma("clang diagnostic push")                                        \
		_Pragma("clang diagnostic ignored \"-Wunused-local-typedefs\"")     \
		_Pragma("clang diagnostic ignored \"-Wunused-variable\"")           \
		_Pragma("clang diagnostic ignored \"-Wunused-parameter\"")          \
		_Pragma("clang diagnostic ignored \"-Wlanguage-extension-token\"")  \
		_Pragma("clang diagnostic ignored \"-Wfloat-equal\"")          

#define CLANG_IGNORE_POP _Pragma("clang diagnostic pop")

 //#ifndef thread_local
 //# if __STDC_VERSION__ >= 201112 && !defined __STDC_NO_THREADS__
 //#  define thread_local _Thread_local
 //# elif defined _WIN32 && ( \
 //	   defined _MSC_VER || \
 //	   defined __ICL || \
 //	   defined __DMC__ || \
 //	   defined __BORLANDC__ )
 //#  define thread_local __declspec(thread) 
 // /* note that ICC (linux) and Clang are covered by __GNUC__ */
 //# elif defined __GNUC__ || \
 //	   defined __SUNPRO_C || \
 //	   defined __xlC__
 //#  define thread_local __thread
 //# else
 //#  error "Cannot define thread_local"
 //# endif
 //#endif

  //#ifdef __STDC_NO_ATOMICS__
  //#error Please use C11 or better with ATOMICS
  //#else
  //#include <stdatomic.h>
  //#endif

#define DBJ_SWAP(x,y) do {\
typeof(x) T_ = x;      \
x = y;                 \
y = T_;                \
 } while (0)

 /* NDEBUG == RELEASE */
#include <assert.h>
/////////////////////////////////////////////////////////////////////////

#undef NOMEM_POLICY

#ifdef NDEBUG 
#define NOMEM_POLICY( BOOLEXP_ ) ((void)BOOLEXP_)
#else // ! NDEBUG == DEBUG
#define NOMEM_POLICY( BOOLEXP_ ) if (! BOOLEXP_ ) { perror( __FILE__ ", Could not allocate memory!"); exit(-1); }
#endif // ! NDEBUG

// for when we are sure ARR is the array
#define DBJ_CNT(ARR) ( sizeof(ARR) / sizeof(ARR[0]) )

#undef MALLOC_WITH_POLICY
#define MALLOC_WITH_POLICY(PTR_ , SIZE_)  do { PTR_ = malloc( SIZE_); NOMEM_POLICY(PTR_); } while(0)

#undef CALLOC_WITH_POLICY
#define CALLOC_WITH_POLICY(PTR_ ,R_,C_, SIZE_)  do { PTR_ = calloc(R_ * C_, SIZE_); NOMEM_POLICY(PTR_); } while(0)

#define DBJ_FREE(P_) do { if (P_){ free(P_); P_ = NULL; }  }while(0)

#undef DBJ_API
#define DBJ_API static

#pragma endregion // common_trash

#pragma region common data 
/////////////////////////////////////////////////////////////////////////
// dimensions
#if DBJ_BENCHMARKING

// NOTE: here we use stack based matrices , thus be carefull with sizes
//       UBENCH repeats execution so matrix size is not the prevailing factor
#define DBJ_MX_A_ROWS 0xF
#define DBJ_MX_A_COLS 0xFF
#define DBJ_MX_B_ROWS DBJ_MX_A_COLS
#define DBJ_MX_B_COLS 0xF

#else // testing
	/*
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

// NOTE: these are compile time typedefs 
// we can create them if we do not use Variably Modified Types (VMT)

typedef dbj_matrix_data_type(*dbj_mx_a_pointer)[DBJ_MX_A_COLS][DBJ_MX_A_ROWS];
typedef dbj_matrix_data_type(*dbj_mx_b_pointer)[DBJ_MX_B_COLS][DBJ_MX_B_ROWS];
typedef dbj_matrix_data_type(*dbj_mx_r_pointer)[DBJ_MX_R_COLS][DBJ_MX_R_ROWS];

typedef dbj_matrix_data_type(*dbj_mx_a_row)[DBJ_MX_A_COLS];
typedef dbj_matrix_data_type(*dbj_mx_b_row)[DBJ_MX_B_COLS];
typedef dbj_matrix_data_type(*dbj_mx_r_row)[DBJ_MX_R_COLS];

#pragma endregion // common data 

#pragma region matrix functions and various matmuls

DBJ_API void* matrix_arr_init
(const unsigned rows_a, const unsigned cols_a, dbj_matrix_data_type a[static rows_a][cols_a]) {
	for (unsigned i = 0; i < rows_a; i++) {
		for (unsigned j = 0; j < cols_a; j++) {
			a[i][j] = (dbj_matrix_data_type)(i * cols_a + j);
		}
	}
	return a;
}

#define dbj_matrix_size_in_bytes( rows_, cols_, type_ ) ( rows_ * cols_ * sizeof(type_) )

DBJ_API void dbj_matrix_transpose(
	const unsigned rows_m,
	const unsigned cols_m,
	const dbj_matrix_data_type m[static rows_m][cols_m],
	dbj_matrix_data_type t[static cols_m][rows_m])
{
	for (size_t i = 0; i < rows_m; i++) {
		for (size_t j = 0; j < cols_m; j++) {
			t[j][i] = m[i][j];
		}
	}
}

DBJ_API dbj_matrix_data_type sdot_1
(int n, const dbj_matrix_data_type x[static n], const dbj_matrix_data_type y[static n])
{
	dbj_matrix_data_type s = (dbj_matrix_data_type)0;
	for (int i = 0; i < n; ++i) s += x[i] * y[i];
	return s;
}

DBJ_API dbj_matrix_data_type sdot_8
(int n, const dbj_matrix_data_type x[static n], const dbj_matrix_data_type y[static n])
{
	int i, n8 = n >> 3 << 3;
	dbj_matrix_data_type s = (dbj_matrix_data_type)0, t[8] = { (dbj_matrix_data_type)0 };
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
	for (s = (dbj_matrix_data_type)0; i < n; ++i) s += x[i] * y[i];
	s += t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
	return s;
}

// the most "by the book" C matrix mutliplication function
// author has added the static keyword for sizes
// this is using VLA/VMT features
// the key fact might be this is the matrix mutliplication so 
// "severley optimized" by compilers there is no point investing
// in finding faster algorithms, including SSE/AVX usage
DBJ_API void the_most_by_the_book_matrix_mult(
	size_t a_rows,
	size_t a_cols,
	size_t b_cols,
	dbj_matrix_data_type A[static a_rows][a_cols],
	dbj_matrix_data_type B[static a_cols][b_cols],
	dbj_matrix_data_type C[static a_rows][b_cols])
{
	for (size_t i = 0; i < a_rows; ++i) {
		for (size_t j = 0; j < b_cols; ++j) {
			C[i][j] = 0.0;
			for (size_t l = 0; l < a_cols; ++l) {
				C[i][j] += A[i][l] * B[l][j];
			}
		}
	}
}


/*
 use 1D aray as matrix type + index calculation of "matrix" [row][col]
 this is in here because it is curiously and persistently the fastest matmul
 */
DBJ_API dbj_matrix_data_type* matmul_mx_as_array
(
	const size_t a_rows, const size_t a_cols, const size_t b_cols,
	dbj_matrix_data_type* a, dbj_matrix_data_type* b, dbj_matrix_data_type* c
)
{
	/*
	the matmul dimensional requirements
	A rows    == B columns
	A columns == B rows
	R rows    == A rows
	R columns == B columns
	*/
	for (size_t i = 0; i < a_rows; i++) {
		for (size_t k = 0; k < b_cols; k++) {
			dbj_matrix_data_type sum = (dbj_matrix_data_type)0.0;
			for (size_t j = 0; j < a_cols /* same as b rows */; j++) {
				sum += a[i * a_cols + j] * b[j * a_rows + k];
			}
			c[i * a_rows + k] = sum;
		}
	}
	return c;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant transposes matrix b, and it's a lot faster. */
DBJ_API dbj_matrix_data_type* matmul_mx_as_array_another
(const size_t a_rows, const size_t a_cols, const size_t b_cols,
	dbj_matrix_data_type* a, dbj_matrix_data_type* b, dbj_matrix_data_type* c
)
{
	// VLA may not be initialized
	dbj_matrix_data_type bt[a_rows * b_cols]; // = { (dbj_matrix_data_type)0 };

	dbj_matrix_transpose(a_rows, b_cols, b, bt);

	for (unsigned i = 0; i < a_rows; i++) {
		for (unsigned k = 0; k < b_cols; k++) {
			dbj_matrix_data_type sum = 0.0;
			for (unsigned j = 0; j < a_cols; j++) {
				sum += a[i * a_cols + j] * bt[k * b_cols + j];
			}
			c[i * b_cols + k] = sum;
		}
	}
	return c;
}

// this is VMT based
DBJ_API void* matmul_transpose_sdot(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	dbj_matrix_data_type a[static n_a_rows][n_a_cols],
	dbj_matrix_data_type b[static n_a_rows][n_b_cols],
	dbj_matrix_data_type m[static n_a_rows][n_b_cols]
)
{
	int  n_b_rows = n_a_cols;

	// variable-sized object may not be initialized
	dbj_matrix_data_type bt[n_a_rows][n_b_cols];

	dbj_matrix_transpose(n_a_rows, n_b_cols, b, bt);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_8(n_a_cols, a[i], bt[j]);
	return m;
}

DBJ_API void* matmul_transpose_sdot_faster(
	const unsigned n_a_rows, const unsigned n_a_cols, const unsigned n_b_cols,
	dbj_matrix_data_type a[static n_a_rows][n_a_cols],
	dbj_matrix_data_type b[static n_a_rows][n_b_cols],
	dbj_matrix_data_type m[static n_a_rows][n_b_cols]
)
{
	int n_b_rows = n_a_cols;

	// variable-sized object may not be initialized
	dbj_matrix_data_type bT[n_a_rows][n_b_cols];
	dbj_matrix_transpose(n_a_rows, n_b_cols, b, bT);

	for (unsigned i = 0; i < n_a_rows; ++i)
		for (unsigned j = 0; j < n_b_cols; ++j)
			m[i][j] = sdot_1(n_a_cols, a[i], bT[j]);
	return m;
}


#pragma endregion // matrix functions and various matmuls

#pragma region common for testing or benchmarking

// ubench functions have no parameters
// thus we use common data aka globals
typedef struct {

	const unsigned rows_a;
	const unsigned cols_a;
	const unsigned rows_b;
	const unsigned cols_b;
	const unsigned rows_r;
	const unsigned cols_r;
	// the matrixes
	dbj_matrix_data_type a[DBJ_MX_A_ROWS][DBJ_MX_A_COLS];
	dbj_matrix_data_type b[DBJ_MX_B_ROWS][DBJ_MX_B_COLS];
	dbj_matrix_data_type r[DBJ_MX_R_ROWS][DBJ_MX_R_COLS]; /* rezult size is a rows * b cols */

} app_data_type;

#define reset_test_result() do { \
dbj_matrix_data_type (*rap)[DBJ_MX_R_ROWS * DBJ_MX_R_COLS] = app_data.r ; \
memset( rap, 0, sizeof(dbj_matrix_data_type[DBJ_MX_R_ROWS * DBJ_MX_R_COLS]));    \
} while (0)

DBJ_API app_data_type app_data = {
			.rows_a = DBJ_MX_A_ROWS,
				.cols_a = DBJ_MX_A_COLS,
				.rows_b = DBJ_MX_B_ROWS,
				.cols_b = DBJ_MX_B_COLS,
				.rows_r = DBJ_MX_A_ROWS, /* the result */
				.cols_r = DBJ_MX_B_COLS,
				// the rest is auto zeroed 
#if !DBJ_BENCHMARKING
// testing data 
			{ {1,2},{3,4} },
			{ {5,6},{7,8} },
			{ {0,0},{0,0} }
#endif // !DBJ_BENCHMARKING
};

DBJ_API void app_start(void)
{
#undef DBJ_APP_KIND

#if DBJ_BENCHMARKING

#define DBJ_APP_KIND  "BENCHMARKING"

	matrix_arr_init(app_data.rows_a, app_data.cols_a, app_data.a);
	matrix_arr_init(app_data.rows_b, app_data.cols_b, app_data.b);
	matrix_arr_init(app_data.rows_r, app_data.cols_r, app_data.r);

#else // TESTING 

#define DBJ_APP_KIND  "TESTING"

	/*
	 *     ! 1 2 |      | 5 6 |       | 19 22 |
	 *     |     |  x   |     |  =    |       |
	 *     | 3 4 |      | 7 8 |       | 43 50 |
	 */
	assert(app_data.rows_a * app_data.cols_a == 4);
	assert(app_data.rows_b * app_data.cols_b == 4);
	assert(app_data.rows_r * app_data.cols_r == 4);

#endif // ! DBJ_BENCHMARKING

	const float size_a = dbj_matrix_size_in_bytes(app_data.rows_a, app_data.cols_a, dbj_matrix_data_type) / 1024.0f;
	const float size_b = dbj_matrix_size_in_bytes(app_data.rows_b, app_data.cols_b, dbj_matrix_data_type) / 1024.0f;
	const float size_r = dbj_matrix_size_in_bytes(app_data.rows_r, app_data.cols_r, dbj_matrix_data_type) / 1024.0f;

	fprintf(stderr, "\n\n" DBJ_VT_RED DBJ_APP_KIND DBJ_VT_RESET " various matrix multiplication algorithms"
		"\n(c) 2021 by dbj dot org, https://dbj.org/license_dbj \nTimestamp: %s"
		"\n\nMatrices are\n"
		"\nA :%4d * %4d * sizeof(%s) == %4.2f KB"
		"\nB :%4d * %4d * sizeof(%s) == %4.2f KB"
		"\nR :%4d * %4d * sizeof(%s) == %4.2f KB\n\n"
		, DBJ_BUILD_TIMESTAMP,
		app_data.rows_a, app_data.cols_a, dbj_matrix_data_type_name, size_a,
		app_data.rows_b, app_data.cols_b, dbj_matrix_data_type_name, size_b,
		app_data.rows_r, app_data.cols_r, dbj_matrix_data_type_name, size_r
	);
#undef DBJ_APP_KIND	
}

DBJ_API void app_end(void)
{
	free(app_data_ptr());
	printf(" " DBJ_VT_RESET " ");

}
/////////////////////////////////////////////////////////////////////////

#if DBJ_BENCHMARKING

// rezult reset and checking are done in UTEST's, see bellow

UBENCH(matmul, matmul_transpose_sdot_faster) {
	matmul_transpose_sdot_faster(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data.a, app_data.b, app_data.r);
}

UBENCH(matmul, matmul_transpose_sdot) {
	matmul_transpose_sdot(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data.a, app_data.b, app_data.r);
}

UBENCH(matmul, matmul_mx_as_array_another) {
	matmul_mx_as_array_another(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
}

UBENCH(matmul, matmul_mx_as_array) {
	matmul_mx_as_array(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
}

UBENCH(matmul, the_most_by_the_book_matrix_mult)
{
	the_most_by_the_book_matrix_mult(
		DBJ_MX_A_ROWS,
		DBJ_MX_A_COLS,
		DBJ_MX_B_COLS,
		app_data.a,
		app_data.b,
		app_data.r
	);
}

/////////////////////////////////////////////////////////////////////////

UBENCH_STATE();

int main(const int argc, const char** argv)
{
#ifdef WIN32
	// VT100 ESC codes kick-start hack
	// 2021-07-10
	// works and necessary 
	// Microsoft Windows [Version 10.0.19042.1052]
	system(" ");
#endif
	app_start();
	return ubench_main(argc, argv);
	app_end();
}

#else // testing /////////////////////////////////////////////////////
/*
 *     ! 1 2 |      | 5 6 |       | 19 22 |
 *     |     |  x   |     |  =    |       |
 *     | 3 4 |      | 7 8 |       | 43 50 |
 */
#define check_test_input() \
do {\
	EXPECT_EQ(app_data.a[0][0] , (dbj_matrix_data_type)1);\
	EXPECT_EQ(app_data.a[0][1] , (dbj_matrix_data_type)2);\
	EXPECT_EQ(app_data.a[1][0] , (dbj_matrix_data_type)3);\
	EXPECT_EQ(app_data.a[1][1] , (dbj_matrix_data_type)4);\
\
	EXPECT_EQ(app_data.b[0][0] , (dbj_matrix_data_type)5);\
	EXPECT_EQ(app_data.b[0][1] , (dbj_matrix_data_type)6);\
	EXPECT_EQ(app_data.b[1][0] , (dbj_matrix_data_type)7);\
	EXPECT_EQ(app_data.b[1][1] , (dbj_matrix_data_type)8);\
} while(0)

#define check_test_result() \
do {\
	EXPECT_EQ(app_data.r[0][0] , (dbj_matrix_data_type)19);\
	EXPECT_EQ(app_data.r[0][1] , (dbj_matrix_data_type)22);\
	EXPECT_EQ(app_data.r[1][0] , (dbj_matrix_data_type)43);\
	EXPECT_EQ(app_data.r[1][1] , (dbj_matrix_data_type)50);\
} while(0)


UTEST(matmul, matmul_transpose_sdot_faster) {
	reset_test_result();
	matmul_transpose_sdot_faster(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data.a, app_data.b, app_data.r);
	check_test_result();
}


UTEST(matmul, matmul_transpose_sdot) {
	reset_test_result();
	matmul_transpose_sdot(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		app_data.a, app_data.b, app_data.r);
	check_test_result();
}

UTEST(matmul, matmul_mx_as_array_another) {
	reset_test_result();
	matmul_mx_as_array_another(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
	check_test_result();
}

UTEST(matmul, matmul_mx_as_array) {
	reset_test_result();
	matmul_mx_as_array(
		DBJ_MX_A_ROWS, DBJ_MX_A_COLS, DBJ_MX_B_COLS,
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
	check_test_result();
}

UTEST(matmul, the_most_by_the_book_matrix_mult) {
	reset_test_result();
	the_most_by_the_book_matrix_mult(
		DBJ_MX_A_ROWS,
		DBJ_MX_B_ROWS,
		DBJ_MX_A_COLS,
		app_data.a,
		app_data.b,
		app_data.r
	);
	check_test_result();
}

UTEST_STATE();

int main(int argc, const char* const argv[]) {
#ifdef WIN32
	// VT100 ESC codes kick-start hack
	// 2021-07-10
	// works and necessary 
	// Microsoft Windows [Version 10.0.19042.1052]
	system(" ");
#endif
	app_start();
	return utest_main(argc, argv);
	app_end();
}

#endif // ! DBJ_BENCHMARKING

#pragma endregion // common for testing or benchmarking
