
#include "../ubench.h/ubench.h"

#include "../dbj_matmul_common.h"

#define DBJ_VARIOUS_SIMPLE_MATMULS_IMPLEMENTATION
#include "simple_matmuls.h"

#undef DBJ_API
#define DBJ_API __attribute__((const)) static


DBJ_API DBJ_MATRIX_DATA_TYPE  get_one(void) { (DBJ_MATRIX_DATA_TYPE)1; }

/**********************************************************************************/

enum {
	max_matrix_side = DBJ_SANITY_MAX_ROW_COUNT,
	testing_matrix_rows = DBJ_MATRIX_SIDE_DIMENSION,
	testing_matrix_cols = DBJ_MATRIX_SIDE_DIMENSION
};

typedef struct {
	unsigned rows;
	unsigned cols;
	DBJ_MATRIX_DATA_TYPE* mx_data_A;
	DBJ_MATRIX_DATA_TYPE* mx_data_B;
	DBJ_MATRIX_DATA_TYPE* mx_data_M;
} test_matmul_common_data;

DBJ_API test_matmul_common_data TMCD_ =
{ .rows = testing_matrix_rows, .cols = testing_matrix_cols, .mx_data_A = NULL, .mx_data_B = NULL, .mx_data_M = NULL };

__attribute__((constructor))
DBJ_API void test_matmul_common_data_constructor(void)
{
	TMCD_.mx_data_A = make_simple_matrix(TMCD_.rows, TMCD_.cols, get_one);
	TMCD_.mx_data_B = make_simple_matrix(TMCD_.rows, TMCD_.cols, get_one);
	TMCD_.mx_data_M = make_simple_matrix(TMCD_.cols, TMCD_.rows, get_one);
}

__attribute__((destructor))
DBJ_API void test_matmul_common_data_destructor(void)
{
	if (TMCD_.mx_data_M) free(TMCD_.mx_data_M);
	if (TMCD_.mx_data_A) free(TMCD_.mx_data_A);
	if (TMCD_.mx_data_B) free(TMCD_.mx_data_B);
}

/*
* calling direct or with table lokup makes no difference in speed
*/
#define UBENCH_COMMON_BODY(DBJ_MATMUL_API_FUN_ID) \
dbj_simple_matmuls_algo_table[DBJ_MATMUL_API_FUN_ID].function(\
	TMCD_.rows, TMCD_.cols, TMCD_.rows, /* BUG?! */ \
	TMCD_.mx_data_A , TMCD_.mx_data_B, TMCD_.mx_data_M \
)

UBENCH(simple, matmul_0) { UBENCH_COMMON_BODY(0); }
UBENCH(simple, matmul_0_0) { 
simple_mat_mul_0_0(
	TMCD_.rows, TMCD_.cols, TMCD_.rows, /* BUG?! */ 
	TMCD_.mx_data_A , TMCD_.mx_data_B, TMCD_.mx_data_M 
);
}

#if 0

UBENCH(simple, matmul_1) { UBENCH_COMMON_BODY(2); }
#if __SSE__
UBENCH(simple, matmul_2_sse) { UBENCH_COMMON_BODY(3); }
UBENCH(simple, matmul_7_sse) { UBENCH_COMMON_BODY(4); }
#endif // __SSE__
UBENCH(simple, matmul_3) { UBENCH_COMMON_BODY(5); }
UBENCH(simple, matmul_4) { UBENCH_COMMON_BODY(6); }

#ifdef HAVE_CBLAS
UBENCH(simple, matmul_5_cblas) { UBENCH_COMMON_BODY(7); }
UBENCH(simple, matmul_6_cblas) { UBENCH_COMMON_BODY(8); }
#endif // HAVE_CBLAS

#endif // 0
