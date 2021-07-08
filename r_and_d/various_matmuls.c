
#include "../ubench.h/ubench.h"

#include "../dbj_matmul_common.h"

//#include <assert.h>
//#include <stdlib.h>
//#include <stdint.h>
//#include <string.h>
//#include <stdio.h>
//#include <time.h>


#define DBJ_FLOAT_VARIOUS_MATMULS_IMPLEMENTATION
#include "dbj_float_various_matmuls.h"
#include "pseudo_random.h"
// #include "dbj_matrix.h"

#undef DBJ_API
#define DBJ_API __attribute__((const)) static

#ifdef NEED_DBJ_SLEEP
/* ----------------------------------------------------------------
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

#endif // NEED_DBJ_SLEEP

/**********************************************************************************/

enum {
	default_multiplication_algorithm = 4, max_matrix_side = DBJ_SANITY_MAX_ROW_COUNT,
	testing_matrix_rows = DBJ_MATRIX_SIDE_DIMENSION,
	testing_matrix_cols = DBJ_MATRIX_SIDE_DIMENSION
};

typedef struct {
	unsigned rows;
	unsigned cols;
	float_matrix_struct* mx_struct_A;
	float_matrix_struct* mx_struct_B;
	float_matrix_struct* mx_struct_M;
} test_matmul_common_data;

DBJ_API test_matmul_common_data TMCD_ = { .rows = testing_matrix_rows, .cols = testing_matrix_cols, .mx_struct_A = NULL, .mx_struct_B = NULL, .mx_struct_M = NULL };

__attribute__((constructor))
DBJ_API void test_matmul_common_data_constructor(void)
{
	TMCD_.mx_struct_A = make_random_float_matrix(TMCD_.rows, TMCD_.cols, pseudo_rand_float);
	TMCD_.mx_struct_B = make_random_float_matrix(TMCD_.rows, TMCD_.cols, pseudo_rand_float);

	TMCD_.mx_struct_M = new_float_matrix(TMCD_.cols, TMCD_.rows);
}

__attribute__((destructor))
DBJ_API void test_matmul_common_data_destructor(void)
{
	if (TMCD_.mx_struct_M) free(TMCD_.mx_struct_M);
	if (TMCD_.mx_struct_A) free(TMCD_.mx_struct_A);
	if (TMCD_.mx_struct_B) free(TMCD_.mx_struct_B);
}

/*
* calling direct or with table lokup makes no difference in speed
*/
#define UBENCH_COMMON_BODY(DBJ_MATMUL_API_FUN_ID) \
dbj_float_matmuls_algo_table[DBJ_MATMUL_API_FUN_ID].function(\
	TMCD_.rows, TMCD_.cols, TMCD_.rows, /* BUG?! */ \
	TMCD_.mx_struct_A->data , TMCD_.mx_struct_B->data, TMCD_.mx_struct_M->data\
)

UBENCH(dbj_various_float_matmuls, matmul_0) { UBENCH_COMMON_BODY(0); }
UBENCH(dbj_various_float_matmuls, matmul_1) { UBENCH_COMMON_BODY(1); }
UBENCH(dbj_various_float_matmuls, matmul_3) { UBENCH_COMMON_BODY(3); }
UBENCH(dbj_various_float_matmuls, matmul_4) { UBENCH_COMMON_BODY(4); }
#ifdef HAVE_CBLAS
UBENCH(dbj_various_float_matmuls, matmul_5_cblas) { UBENCH_COMMON_BODY(5); }
UBENCH(dbj_various_float_matmuls, matmul_6_cblas) { UBENCH_COMMON_BODY(6); }
#endif // HAVE_CBLAS
#if __SSE__
UBENCH(dbj_various_float_matmuls, matmul_2_sse) { UBENCH_COMMON_BODY(2); }
UBENCH(dbj_various_float_matmuls, matmul_7_sse) { UBENCH_COMMON_BODY(7); }
#endif // __SSE__
