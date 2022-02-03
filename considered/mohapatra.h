#pragma once

#define LOOP(C, R) for (unsigned C = 0; C < R; ++C)

/*
https://raw.githubusercontent.com/ShrohanMohapatra/matrix_multiply_quadratic/master/matrix_multiply_test.py

Original Python code by https://github.com/ShrohanMohapatra

Copyright (c) 2019, ShrohanMohapatra
All rights reserved.

transformation to C -- (c) 20201 by dbj@dbj.org
 */

static void mohapatra(const unsigned Arows, const unsigned Acols, /* == Brows */
                      const unsigned Bcols, double A[Arows][Acols],
                      double B[Acols][Bcols], double E[Arows][Bcols]) {
  double maxi = 0;
  LOOP(i, Arows)
  LOOP(j, Bcols) {
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

  LOOP(i, Arows) {
    double sum_1 = 0;
    LOOP(j, Acols)
    sum_1 = sum_1 * (pow(10, P)) + A[i][j];
    C[i] = sum_1;
  }
  LOOP(j, Acols) {
    double sum_1 = 0;
    LOOP(i, Arows) { sum_1 = sum_1 * (pow(10, P)) + B[Arows - 1 - i][j]; }
    D[j] = sum_1;
  }

  LOOP(i, Arows) {
    LOOP(j, Bcols) {
      E[i][j] =
          (int)(C[i] * D[j] / (pow(10, (P * (Arows - 1))))) % (int)(pow(10, P));
    }
  }
} // mohapatra

// check mohapatra one more
UTEST(matmul, mohapatra_separate_test) {
  dbj_matrix_data_type a[6][3] = {{1, 2, 3}, {1, 2, 3}, {1, 2, 3},
                                  {1, 2, 3}, {1, 2, 3}, {1, 2, 3}};
  dbj_matrix_data_type b[3][6] = {
      {1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
  dbj_matrix_data_type r1[6][6] = {0};
  dbj_matrix_data_type r2[6][6] = {0};

  mohapatra(6, 3, a, b, r1);
  ijk_matmul(6, 3, 6, a, b, r2);
  EXPECT_TRUE(compare_matrices(6, 6, r1, r2));

  {
    printf("\nR1");
    print_matrix("%4.2f", 6, 6, r1);
    printf("\nR2");
    print_matrix("%4.2f", 6, 6, r2);
    printf("\n");
  }
}

UTEST(matmul, mohapatra) {
  reset_test_result(app_data);
  mohapatra(DBJ_MX_A_ROWS, DBJ_MX_A_COLS, app_data->a, app_data->b,
            app_data->r);
  check_test_result();
}