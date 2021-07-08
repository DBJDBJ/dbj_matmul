#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <time.h>


#define DBJ_FLOAT_VARIOUS_MATMULS_IMPLEMENTATION
#include "dbj_float_various_matmuls.h"
#include "pseudo_random.h"
// #include "dbj_matrix.h"

#undef DBJ_API
#define DBJ_API static

#pragma region dbj_sleep
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

#pragma endregion dbj_sleep

/**********************************************************************************/

enum {
	default_multiplication_algorithm = 4, max_matrix_side = DBJ_SANITY_MAX_ROW_COUNT,
	testing_matrix_rows = 0xFF,
	testing_matrix_cols = 0xFF
};

typedef struct {
	unsigned rows;
	unsigned cols;
	float_matrix_struct* mx_struct_A;
	float_matrix_struct* mx_struct_B;
	float_matrix_struct* mx_struct_M;
} test_matmul_common_data;

DBJ_API test_matmul_common_data TMCD_ = { .rows = testing_matrix_rows, .cols = testing_matrix_cols, NULL, NULL, NULL };

__attribute__((constructor))
DBJ_API test_matmul_common_data_constructor(void)
{
	TMCD_.mx_struct_A = make_random_float_matrix(TMCD_.rows, TMCD_.cols, pseudo_rand_float);
	TMCD_.mx_struct_B = make_random_float_matrix(TMCD_.rows, TMCD_.cols, pseudo_rand_float);

	TMCD_.mx_struct_M = new_float_matrix(TMCD_.cols, TMCD_.rows);
}

__attribute__((destructor))
DBJ_API test_matmul_common_data_destructor(void)
{
	if (TMCD_.mx_struct_M) free(TMCD_.mx_struct_M);
	if (TMCD_.mx_struct_A) free(TMCD_.mx_struct_A);
	if (TMCD_.mx_struct_B) free(TMCD_.mx_struct_B);
}


/*
These are square matrices, matrix_side_length is the length of the side of that square
matrix data type is pointer to pointer
matrix element type is float
*/
DBJ_API int test_matmul(unsigned algorithm_id, description_function_pair dfp)
{
	int algo = algorithm_id; (void)algo;
	const char* algo_name = dfp.description;
	mat_mul_function mm_fun = dfp.function;

	assert(algo_name);

	DBJ_MATRIX_ALIAS(matrix, float, TMCD_.cols);

	//float_matrix_struct* mx_struct_A = NULL, * mx_struct_B = NULL, * mx_struct_M = NULL;

	//assert(matrix_side_length < max_matrix_side);
	//// release mode: quietly adjusting to max side .. not exactly good but ok in this context
	//matrix_side_length = matrix_side_length % (size_t)max_matrix_side;

	//mx_struct_A = make_random_float_matrix(matrix_side_length, matrix_side_length, pseudo_rand_float);
	//mx_struct_B = make_random_float_matrix(matrix_side_length, matrix_side_length, pseudo_rand_float);

	//mx_struct_M = new_float_matrix(matrix_side_length, matrix_side_length);

	DBJ_MATRIX_CAST(mx_A, matrix, TMCD_.mx_struct_A);
	DBJ_MATRIX_CAST(mx_B, matrix, TMCD_.mx_struct_B);
	DBJ_MATRIX_CAST(mx_M, matrix, TMCD_.mx_struct_M);


	clock_t start_time_ = clock();

	(void)mm_fun(
		TMCD_.rows, TMCD_.cols, TMCD_.rows, /* BUG?! */
		mx_A, mx_B, mx_M
	);

	fprintf(stderr, "\nAlgorithm %-36s: ", algo_name);
	fprintf(stderr, "\nCPU time: %2.3g sec", (double)(clock() - start_time_) / CLOCKS_PER_SEC);
	// fprintf(stderr, "\nCentral cell: %g", mx_M[matrix_side_length / 2][matrix_side_length / 2]);

	//if (mx_struct_M) free(mx_struct_M);
	//if (mx_struct_A) free(mx_struct_A);
	//if (mx_struct_B) free(mx_struct_B);

	return EXIT_SUCCESS;
}

static void test_various_matmuls(void)
{
	fprintf(stderr, "\n------------------------------------------------");
	fprintf(stderr, "\nMatrix width == height == size: %4d", TMCD_.rows);

	static const unsigned algo_table_count = sizeof(algo_table) / sizeof(algo_table[0]);

	for (unsigned k = 0; k < algo_table_count; ++k)
	{
		fprintf(stderr, "\n------------------------------------------------");
		test_matmul(k, algo_table[k]);

		dbj_sleep(dbj_one_sec);
	}

}

int various_matmuls(const int argc, const char** argv)
{
	(void)argc;
	(void)argv;
	//float_matrix_struct* fmtp = make_random_float_matrix(3, 4);
	//DBJ_MATRIX_FREE(fmtp);
	test_various_matmuls();
	return 42;
}