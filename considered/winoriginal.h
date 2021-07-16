#ifndef WINORIGINAL_H
#define WINORIGINAL_H
//
//  winograd.c
//  strassen-winograd
//
//  Created by Pavel Kravets on 04.10.13.
//  Copyright (c) 2013 Pavel Kravets. All rights reserved.
//

#ifndef WINORIGINAL_SQUARE_SIDE
#error WINORIGINAL_SQUARE_SIDE not defined?
#endif  // WINORIGINAL_SQUARE_SIDE

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
    unsigned colNum;
    unsigned rowNum;
    double data[WINORIGINAL_SQUARE_SIDE * 2] ;
} matrix ;

typedef double(*winoriginal_data_row_ptr)[WINORIGINAL_SQUARE_SIDE] ;

// VMT
#define element(m,x,y) ((winoriginal_data_row_ptr)m->data)[x][y] 

 // set_element(res, i, j, -row[i]-col[j]);
 #define set_element(M_, i, j, V_) element(M_,i,j) = V_ 

static inline void original_winograd_preprocess(matrix* m1, matrix* m2, double* row, double* col)
{
    int d = m1->colNum/2;
    
    for (int i=0; i<m1->rowNum; i++)
    {
        row[i] = 0;
        for (int j=0; j<d; j++)
        {
            row[i] += element(m1, i, 2*j) * element(m1, i, 2*j+1);
        }
    }
    
    for (int i=0; i<m2->colNum; i++)
    {
        col[i] = 0;
        for (int j=0; j<d; j++)
        {
            col[i] += element(m2, 2*j, i) * element(m2, 2*j+1, i);
        }
    }
    
}

static inline void original_winograd_mult(matrix* m1, matrix* m2, matrix* res)
{
    assert (m1->colNum == m2->rowNum);
    assert( res->rowNum == m1->rowNum);
    assert(res->colNum == m2->colNum);
    
    double* row = calloc(m1->rowNum , sizeof(double));
    double* col = calloc(m2->colNum , sizeof(double));
    
    original_winograd_preprocess(m1, m2, row, col);

    int d = m1->colNum/2;
    
    for (int i=0; i<res->rowNum; i++)
    {
        for (int j=0; j<res->colNum; j++)
        {
            set_element(res, i, j, -row[i]-col[j]);
            for (int k=0; k<d; k++)
            {
                set_element(res, i, j, (element(m1, i, 2*k) + element(m2, 2*k+1, j)) * (element(m1, i, 2*k+1)+element(m2, 2*k, j)));
            }
            if (m1->colNum % 2 !=0)
            {
                set_element(res, i, j, element(m1, i, m1->colNum-1)*element(m2, m2->rowNum-1, j));
            }
        }
    }

    free(row);
    free(col);
    
}

#define WINORIGINAL_TEST_OR_BENCH
#ifdef WINORIGINAL_TEST_OR_BENCH

	typedef struct winoriginal_struct {
		unsigned square_side;
		matrix m1; 
		matrix m2; 
		matrix res; 
	} winoriginal_app_data_type ;

    static winoriginal_app_data_type * winoriginal_app_data = 0 ;

__attribute__((constructor))
    static inline void winoriginal_start (void )
{
    winoriginal_app_data = calloc(1, sizeof(winoriginal_app_data_type));    assert(winoriginal_app_data) ;

    winoriginal_app_data->m1.colNum = winoriginal_app_data->m1.rowNum = 
    winoriginal_app_data->m2.colNum = winoriginal_app_data->m2.rowNum = 
    winoriginal_app_data->res.colNum = winoriginal_app_data->res.rowNum = WINORIGINAL_SQUARE_SIDE;

	/*
		 In case of testing we use this constelation
		 of matrices, to check the correctness of algorithms

	 *     ! 1 2 |      | 5 6 |       | 19 22 |
	 *     |     |  x   |     |  =    |       |
	 *     | 3 4 |      | 7 8 |       | 43 50 |
	 */

    winoriginal_app_data->square_side = WINORIGINAL_SQUARE_SIDE;

    winoriginal_app_data->m1.data[0] = (dbj_matrix_data_type)1;
    winoriginal_app_data->m1.data[1] = (dbj_matrix_data_type)2;
    winoriginal_app_data->m1.data[2] = (dbj_matrix_data_type)3;
    winoriginal_app_data->m1.data[3] = (dbj_matrix_data_type)4;

    winoriginal_app_data->m2.data[0] = (dbj_matrix_data_type)5;
    winoriginal_app_data->m2.data[1] = (dbj_matrix_data_type)6;
    winoriginal_app_data->m2.data[2] = (dbj_matrix_data_type)7;
    winoriginal_app_data->m2.data[3] = (dbj_matrix_data_type)8;
}

__attribute__((destructor)) static inline void winoriginal_end ( void ){ 
    free(winoriginal_app_data);
}

#define winoriginal_check_test_result() \
do {\
	EXPECT_EQ(winoriginal_app_data->res.data[0] , (dbj_matrix_data_type)19);\
	EXPECT_EQ(winoriginal_app_data->res.data[1] , (dbj_matrix_data_type)22);\
	EXPECT_EQ(winoriginal_app_data->res.data[2] , (dbj_matrix_data_type)43);\
	EXPECT_EQ(winoriginal_app_data->res.data[3] , (dbj_matrix_data_type)50);\
} while(0)

#define winoriginal_reset_test_result() do { \
memset( winoriginal_app_data->res.data, 0, sizeof(dbj_matrix_data_type[WINORIGINAL_SQUARE_SIDE * WINORIGINAL_SQUARE_SIDE]));    \
} while (0)

UTEST( winograd_original, testing )
{
    winoriginal_reset_test_result();
    original_winograd_mult(
        & winoriginal_app_data->m1,  & winoriginal_app_data->m2,  & winoriginal_app_data->res 
    );
    winoriginal_check_test_result();
}

#undef winoriginal_check_test_result
#undef winoriginal_reset_test_result

#endif // WINORIGINAL_TEST_OR_BENCH

#undef element
#undef set_element
#undef add_to_element_2
#undef add_to_element_4

#endif // WINORIGINAL_H
