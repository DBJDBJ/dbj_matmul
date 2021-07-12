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

#pragma region common_trash

#define CLANG_IGNORE_PUSH \
	_Pragma("clang diagnostic push")                                           \
		_Pragma("clang diagnostic ignored \"-Wunused-local-typedefs\"")     \
		_Pragma("clang diagnostic ignored \"-Wunused-variable\"")     \
		_Pragma("clang diagnostic ignored \"-Wunused-parameter\"")     \
		_Pragma("clang diagnostic ignored \"-Wlanguage-extension-token\"")     \
		_Pragma("clang diagnostic ignored \"-Wfloat-equal\"")          

#define CLANG_IGNORE_POP _Pragma("clang diagnostic pop")

#ifndef thread_local
# if __STDC_VERSION__ >= 201112 && !defined __STDC_NO_THREADS__
#  define thread_local _Thread_local
# elif defined _WIN32 && ( \
	   defined _MSC_VER || \
	   defined __ICL || \
	   defined __DMC__ || \
	   defined __BORLANDC__ )
#  define thread_local __declspec(thread) 
 /* note that ICC (linux) and Clang are covered by __GNUC__ */
# elif defined __GNUC__ || \
	   defined __SUNPRO_C || \
	   defined __xlC__
#  define thread_local __thread
# else
#  error "Cannot define thread_local"
# endif
#endif

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

#if ! DBJ_ON_GODBOLT
#include "build_time_stamp.inc" // DBJ_BUILD_TIMESTAMP 
#   if DBJ_BENCHMARKING
#include "ubench.h/ubench.h"
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
(dbj_matrix_data_type* a, const unsigned rows_a, const unsigned cols_a) {
	for (unsigned i = 0; i < rows_a; i++) {
		for (unsigned j = 0; j < cols_a; j++) {
			a[i * cols_a + j] = (dbj_matrix_data_type)(i * cols_a + j);
		}
	}
	return a;
}

DBJ_API float dbj_matrix_size_in_bytes(const unsigned rows_, const unsigned cols_, size_t data_size_)
{
	return (rows_ * cols_ * (unsigned)data_size_);
}

/*  Takes and returns a new matrix, t, which is a dbj_matrix_transpose of the original one, m.
   It's also flat in memory, i.e., 1-D, but it should be looked at as a matrix
   of m, meaning, rows_t == cols_m, and cols_t == rows_m.
   The original matrix m stays intact.

   NOTE: VMT make this much safr, but no VMT here
*/
DBJ_API void* dbj_matrix_transpose(
	const unsigned rows_m, const unsigned cols_m,
	const dbj_matrix_data_type* m,
	dbj_matrix_data_type* t)
{
	for (size_t i = 0; i < rows_m; i++) {
		for (size_t j = 0; j < cols_m; j++) {
			t[j * rows_m + i] = m[i * cols_m + j];
		}
	}
	return t;
}

// text bool matmul with VMT arguments
// rezult matrix m is sent in as pre allocated
// notice the dimensions requirements for a,b, and m
DBJ_API void* matmul_row_pointers
(
	dbj_mx_a_row a, dbj_mx_b_row b, dbj_mx_r_row m
)
{
	for (unsigned c = 0; c < DBJ_MX_A_ROWS; c++) {
		for (unsigned d = 0; d < DBJ_MX_B_COLS; d++) {
			for (unsigned k = 0; k < DBJ_MX_B_ROWS; k++) {
				m[c][d] += a[c][k] * b[k][d];
			}
		}
	}
	return m;
}

/*
 use 1D aray as matrix type + index calculation of "matrix" [row][col]
 */
DBJ_API dbj_matrix_data_type* matmul_mx_as_array
(dbj_matrix_data_type* a, dbj_matrix_data_type* b, dbj_matrix_data_type* c)
{

	for (size_t i = 0; i < DBJ_MX_A_ROWS; i++) {
		for (size_t k = 0; k < DBJ_MX_B_COLS; k++) {
			dbj_matrix_data_type sum = (dbj_matrix_data_type)0.0;
			for (size_t j = 0; j < DBJ_MX_B_ROWS; j++) {
				sum += a[i * DBJ_MX_A_COLS + j] * b[j * DBJ_MX_B_COLS + k];
			}
			c[i * DBJ_MX_B_COLS + k] = sum;
		}
	}

	return c;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant transposes matrix b, and it's a lot faster. */
DBJ_API
dbj_matrix_data_type* matmul_mx_as_array_faster
(dbj_matrix_data_type* a, dbj_matrix_data_type* b, dbj_matrix_data_type* c)
{

	dbj_matrix_data_type bt[DBJ_MX_A_ROWS * DBJ_MX_B_COLS] = { (dbj_matrix_data_type)0 };

	dbj_matrix_transpose(DBJ_MX_A_ROWS, DBJ_MX_B_COLS, b, bt);

	for (unsigned i = 0; i < DBJ_MX_A_ROWS; i++) {
		for (unsigned k = 0; k < DBJ_MX_B_COLS; k++) {
			dbj_matrix_data_type sum = 0.0;
			for (unsigned j = 0; j < DBJ_MX_A_COLS; j++) {
				sum += a[i * DBJ_MX_A_COLS + j] * bt[k * DBJ_MX_B_COLS + j];
			}
			c[i * DBJ_MX_B_COLS + k] = sum;
		}
	}
	return c;
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
	dbj_matrix_data_type a[DBJ_MX_A_ROWS * DBJ_MX_A_COLS];
	dbj_matrix_data_type b[DBJ_MX_B_ROWS * DBJ_MX_B_COLS];
	dbj_matrix_data_type r[DBJ_MX_R_ROWS * DBJ_MX_R_COLS]; /* rezult size is a rows * b cols */

} app_data_type;

DBJ_API app_data_type app_data = {
	.rows_a = DBJ_MX_A_ROWS ,
	.cols_a = DBJ_MX_A_COLS ,
	.rows_b = DBJ_MX_B_ROWS ,
	.cols_b = DBJ_MX_B_COLS ,
	.rows_r = DBJ_MX_A_ROWS , /* the result */
	.cols_r = DBJ_MX_B_COLS ,
	// the rest is auto zeroed 
	#if !DBJ_BENCHMARKING
	// testing data 
	{ 1,2,3,4 },
	{ 5,6,7,8 },
	{ 0,0,0,0 }
	#endif // !DBJ_BENCHMARKING
};

DBJ_API void app_start(void)
{
#if DBJ_BENCHMARKING

#define DBJ_APP_KIND  "BENCHMARKING"

	matrix_arr_init(app_data.a, app_data.rows_a, app_data.cols_a);
	matrix_arr_init(app_data.b, app_data.rows_b, app_data.cols_b);
	matrix_arr_init(app_data.r, app_data.rows_a, app_data.cols_b);

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

	fprintf(stderr, "\n\n" DBJ_APP_KIND " various matrix multiplication algorithms"
		"\n(c) 2021 by dbj dot org, https://dbj.org/license_dbj , timestamp: %s"
		"\nMatrices are\n"
		"\nA : %d * %d * sizeof(%s) == %.2f KB"
		"\nB : %d * %d * sizeof(%s) == %.2f KB"
		"\nR : %d * %d * sizeof(%s) == %.2f KB\n\n"
		, DBJ_BUILD_TIMESTAMP,
		app_data.rows_a, app_data.cols_a, dbj_matrix_data_type_name, dbj_matrix_size_in_bytes(app_data.rows_a, app_data.cols_a, sizeof(dbj_matrix_data_type)) / 1024.0f,
		app_data.rows_b, app_data.cols_b, dbj_matrix_data_type_name, dbj_matrix_size_in_bytes(app_data.rows_b, app_data.cols_b, sizeof(dbj_matrix_data_type)) / 1024.0f,
		app_data.rows_r, app_data.cols_r, dbj_matrix_data_type_name, dbj_matrix_size_in_bytes(app_data.rows_r, app_data.cols_r, sizeof(dbj_matrix_data_type)) / 1024.0f
	);
#undef DBJ_APP_KIND	
}

DBJ_API void app_end(void)
{

}
/////////////////////////////////////////////////////////////////////////

#if DBJ_BENCHMARKING

// rezult reset and checking are done in UTEST's, see bellow

UBENCH(matmul, matmul_mx_as_array_faster) {
	matmul_mx_as_array_faster(
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
}

UBENCH(matmul, matmul_mx_as_array) {
	matmul_mx_as_array(
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
}

UBENCH(matmul, matmul_row_ptr) {
	matmul_row_pointers(
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
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
	EXPECT_EQ(app_data.a[0] , (dbj_matrix_data_type)1);\
	EXPECT_EQ(app_data.a[1] , (dbj_matrix_data_type)2);\
	EXPECT_EQ(app_data.a[2] , (dbj_matrix_data_type)3);\
	EXPECT_EQ(app_data.a[3] , (dbj_matrix_data_type)4);\
\
	EXPECT_EQ(app_data.b[0] , (dbj_matrix_data_type)5);\
	EXPECT_EQ(app_data.b[1] , (dbj_matrix_data_type)6);\
	EXPECT_EQ(app_data.b[2] , (dbj_matrix_data_type)7);\
	EXPECT_EQ(app_data.b[3] , (dbj_matrix_data_type)8);\
} while(0)

#define check_test_result() \
do {\
	EXPECT_EQ(app_data.r[0] , (dbj_matrix_data_type)19);\
	EXPECT_EQ(app_data.r[1] , (dbj_matrix_data_type)22);\
	EXPECT_EQ(app_data.r[2] , (dbj_matrix_data_type)43);\
	EXPECT_EQ(app_data.r[3] , (dbj_matrix_data_type)50);\
} while(0)


#define reset_test_result() do { \
for (unsigned k = 0; k < (DBJ_MX_R_COLS * DBJ_MX_R_ROWS); ++k) app_data.r[k] = (dbj_matrix_data_type)0; \
} while (0)

UTEST(matmul, matmul_mx_as_array_faster) {
	reset_test_result();
	matmul_mx_as_array_faster(
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
	check_test_result();
}

UTEST(matmul, matmul_mx_as_array) {
	reset_test_result();
	matmul_mx_as_array(
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
	);
	check_test_result();
}

UTEST(matmul, matmul_row_ptr) {
	reset_test_result();
	matmul_row_pointers(
		(void*)app_data.a,
		(void*)app_data.b,
		(void*)app_data.r
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
