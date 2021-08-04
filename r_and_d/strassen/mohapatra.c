
#include "dbj_matrix.h" // Matrix

// https://godbolt.org/z/7zsPEG6EE
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*
         In case of testing we use this constelation
         of matrices, to check the correctness of algorithms

 *     ! 1 2 |      | 5 6 |       | 19 22 |
 *     |     |  x   |     |  =    |       |
 *     | 3 4 |      | 7 8 |       | 43 50 |
 */
enum { N = 2 };

#ifndef FOR
#define FOR(C, R) for (unsigned C = 0; C < R; ++C)
#endif

/*
Original algorithm:

https://raw.githubusercontent.com/ShrohanMohapatra/matrix_multiply_quadratic/master/matrix_multiply_test.py

BSD 3-Clause License

Copyright (c) 2019, ShrohanMohapatra
All rights reserved.

Transformation to C (c) 2021 by dbj@dbj.org
 */

static void mohapatra(double A[N][N], double B[N][N], double E[N][N]) {
  double maxi = 0;
  FOR(i, N)
  FOR(j, N) {
    if (maxi < A[i][j])
      maxi = A[i][j];
    if (maxi < B[i][j])
      maxi = B[i][j];
  }
  int M = (int)(log10(maxi)) + 1;
  int P = (int)(log10((pow(10, (2 * M)) - 1) * N)) + 1;

  double C[N] = {0}, D[N] = {0};

  FOR(i, N) {
    double sum_1 = 0;
    FOR(j, N) sum_1 = sum_1 * (pow(10, P)) + A[i][j];
    C[i] = sum_1;
  }
  FOR(j, N) {
    double sum_1 = 0;
    FOR(i, N) { sum_1 = sum_1 * (pow(10, P)) + B[N - 1 - i][j]; }
    D[j] = sum_1;
  }

  FOR(i, N) {
    FOR(j, N) {
      E[i][j] =
          (int)(C[i] * D[j] / (pow(10, (P * (N - 1))))) % (int)(pow(10, P));
    }
  }
} // mohapatra

static void printmat(double mx[N][N]) {
  FOR(j, N) {
    FOR(k, N) { printf(" %4.2f ", mx[j][k]); }
    printf("\n");
  }
}

#define MATRIX_A_SIDE 2
#define MATRIX_B_SIDE MATRIX_A_SIDE

static inline Matrix *matrix_setter(Matrix *, unsigned, unsigned);

static inline Matrix *matrix_setter(Matrix *mx, unsigned R, unsigned C) {
  MXY(mx, R, C) = 42;
  return mx;
}

static inline void mohapatra_matrix_using() {
  Matrix *a = matrix_new(MATRIX_A_SIDE, MATRIX_B_SIDE);
  Matrix *b = matrix_new(MATRIX_A_SIDE, MATRIX_B_SIDE);
  Matrix *r = matrix_new(MATRIX_A_SIDE, MATRIX_B_SIDE);

  matrix_foreach(a, matrix_setter);
  matrix_foreach(b, matrix_setter);
  matrix_foreach(r, matrix_setter);

  mohapatra((void *)a->values, (void *)b->values, (void *)r->values);
  printf("\nAfter R = mohapatra(A,B) the result is:\n\n");
  printmat((void *)r->values);
}

/////////////////////////////////////////////////////////////////////////////////////////
#include <fenv.h>
#include <corecrt_math.h>
#include <math.h>
#include <errno.h>
#pragma STDC FENV_ACCESS ON
static void check_matherr(void) {
  printf("MATH_ERRNO is %s\n",
         math_errhandling & MATH_ERRNO ? "set" : "not set");
  printf("MATH_ERREXCEPT is %s\n",
         math_errhandling & MATH_ERREXCEPT ? "set" : "not set");
  feclearexcept(FE_ALL_EXCEPT);
  errno = 0;
  printf("log(0) = %f\n", log(0));
  if (errno == ERANGE)
    perror("errno == ERANGE");
  if (fetestexcept(FE_DIVBYZERO))
    puts("FE_DIVBYZERO (pole error) reported");
}
/////////////////////////////////////////////////////////////////////////////////////////

int main(void) {

  check_matherr();

  mohapatra_matrix_using();

#if 0
  double a[2][2] = {{1, 2}, {3, 4}};
  double b[2][2] = {{5, 6}, {7, 8}};

  //   double row[2] = {0};
  //   double col[2] = {0};

  printf("\nThe constelation"
         "\nof matrices, to check the correctness of algorithms:"
         "\n"
         "\n     | 1 2 |      | 5 6 |       | 19 22 |"
         "\n     |     |  x   |     |  =    |       |"
         "\n     | 3 4 |      | 7 8 |       | 43 50 |"
         "\n");

  {
    double r[2][2] = {{0, 0}, {0, 0}};
    mohapatra(a, b, r);
    printf("\nAfter R = mohapatra(A,B) the result is:\n\n");
    printmat(r);
  }

  printf("\n\nOriginal algorithm:"
         "\n"
         "\nhttps://raw.githubusercontent.com/ShrohanMohapatra/"
         "matrix_multiply_quadratic/master/matrix_multiply_test.py"
         "\n"
         "\nBSD 3-Clause License"
         "\nCopyright (c) 2019, ShrohanMohapatra"
         "\nAll rights reserved."
         "\n"
         "\nTransformation to C (c) 2021 by dbj@dbj.org");
#endif // 0
  return 42;
}
