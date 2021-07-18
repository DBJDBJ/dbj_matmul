
#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

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

#if (defined(__clang__) || defined(__GNUC__))
#define DBJ_CLANGNUC 1
#else
#define DBJ_CLANGNUC 0
#endif

#if !DBJ_ON_GODBOLT
#include "build_time_stamp.inc" // DBJ_BUILD_TIMESTAMP
#if DBJ_BENCHMARKING
#include "ubench/ubench.h"
#else
#include "utest/utest.h"
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
#endif

/*
too much for WIN10 8GB RAM
#define DATA_SIZE 0xFFFFFFFF
*/
#define DATA_SIZE 0xFFFFFFF
#define N_THREADS 2
#define range DATA_SIZE / N_THREADS

static struct {
  double *X; // [DATA_SIZE];
  double *Y; // [DATA_SIZE];
} *app_data = 0;

static void app_start(int argc, const char *const argv[]) {
  app_data = calloc(1, sizeof(*app_data));
  assert(app_data);

  app_data->X = calloc(DATA_SIZE, sizeof(double));
  assert(app_data->X);

  app_data->Y = calloc(DATA_SIZE, sizeof(double));
  assert(app_data->Y);

  FOR(i, DATA_SIZE) {
    app_data->X[i] = 1.0;
    app_data->Y[i] = 2.0;
  }

  printf("\nStarting %s", argv[0]);
  printf("\nMeasuring: double[%d] /= double[%d]\n\n", DATA_SIZE, DATA_SIZE);
}

static void app_end(void) {
  free(app_data->X);
  free(app_data->Y);
  free(app_data);
}

#ifdef _MSC_VER
#pragma range benches
#endif

UBENCH(omp_measuring, adding_two_arrays_of_doubles) {
#pragma omp parallel for num_threads(N_THREADS)
  FOR(i, DATA_SIZE) { app_data->X[i] /= app_data->Y[i]; }
}

UBENCH(NOT_omp_measuring, adding_two_arrays_of_doubles) {
  FOR(i, DATA_SIZE) { app_data->X[i] /= app_data->Y[i]; }
}

#ifdef _MSC_VER
#pragma endrange // benches
#endif

UBENCH_STATE();

int main(int argc, const char *const argv[]) {
  app_start(argc, argv);
  return ubench_main(argc, argv);
  app_end();
}