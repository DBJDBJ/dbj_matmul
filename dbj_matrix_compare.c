#define TESTING 1
#define SAME_ARRAYS 1

#if TESTING
// #include "https://raw.githubusercontent.com/sheredom/utest.h/master/utest.h"
#include "utest/utest.h"
#define BEFORE_MAIN UTEST_STATE
#define MAIN_FUNCTION utest_main
#define BENCH_OR_TEST UTEST
#else
#include "ubench/ubench.h"
#define BEFORE_MAIN UBENCH_STATE
#define MAIN_FUNCTION ubench_main
#define BENCH_OR_TEST UBENCH
#endif // ! TESTING

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// for simple for's
#define FOR(C, N) for (unsigned C = 0; C < (unsigned)N; ++C)

#define DATA_SIZE 100000

typedef double data_type;

static struct {

#if TESTING
  data_type x_d[6];
  data_type y_d[6];
#else
  data_type x_d[DATA_SIZE];
  data_type y_d[DATA_SIZE];
#endif

} *app_data = 0;

static void app_start(int argc, const char *const argv[]) {
  (void)argc;
  (void)argv;
  app_data = calloc(1, sizeof(*app_data));
  assert(app_data);

  FOR(i, DATA_SIZE) {
    app_data->x_d[i] = 1.0;
#if SAME_ARRAYS
    app_data->y_d[i] = 1.0;
#else
    app_data->y_d[i] = 2.0;
#endif // SAME_ARRAYS
  }

  printf("\n\nStarting %s\n", argv[0]);
#if SAME_ARRAYS
#define constelation "SAME"
#else
#define constelation "DIFFERENT"
#endif // SAME_ARRAYS
  printf("\n2 arrays prepared are %s\n\n", constelation);

#undef constelation
}

////////////////////////////////////////////////////////////////

#define close_enough_(a, b)                                                    \
  (fabs((a) - (b)) < 1e-10 * (fabs(a) + fabs(b))) ? 1 : 0

static inline int compare_arrays(const int N, data_type a[static N],
                                 data_type b[static N]) {
  for(unsigned k = 0; k < N; ++k) {
    if (!close_enough_(a[k], b[k]))
      return 0;
  }
  return 1;
}

static inline int compare_matrices(const int N, const int M,
                                   data_type a[static N][M],
                                   data_type b[static N][M]) {
  return compare_arrays(N * M, (data_type *)a, (data_type *)b);
}

#undef close_enough_

////////////////////////////////////////////////////////////////

#if TESTING
UTEST(testing, comparing_two_arrays) {
  EXPECT_TRUE(
      compare_matrices(3, 2, (void *)app_data->x_d, (void *)app_data->y_d));
}
#else
UBENCH(measuring, comparing_two_arrays) {
  (void)compare_arrays(DATA_SIZE, app_data->x_d, app_data->y_d);
}
#endif // ! TESTING

static void app_end(void) { free(app_data); }

BEFORE_MAIN();

int main(int argc, const char *const argv[]) {
#if defined(_WIN32)
  system(" ");
#endif
  app_start(argc, argv);
  return MAIN_FUNCTION(argc, argv);
  app_end();
}