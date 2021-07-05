
#define DBJ_MATMUL_IMPLEMENTATION
#include "dbj_matmul.h"

#if 0

#ifdef __clang__
#include "macro.h"
#endif

#define MATMUL_OUTSTREAM stderr
#define MATMUL_PRINT(...) fprintf( MATMUL_OUTSTREAM, __VA_ARGS__ )

// works server side too
#define NOMEM_POLICY( BOOLEXP_ )\
    if (! BOOLEXP_ ) {\
        perror( __FILE__ ", Could not allocate memory!");\
        exit(-1);\
    }

#define ALLOC_WITH_POLICY(PTR_ , SIZE_)    \
do { \
PTR_ = calloc(1,SIZE_);\
NOMEM_POLICY(PTR_) ;\
} while(0)

/* Initializes vector or matrix, sequentially, with indices. */
_const_ static
void init_seq(double *, const unsigned , const unsigned );
/* Initializes vector or matrix, randomly. */
_const_  static
void init_rand(double *, const unsigned , const unsigned );
/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
    It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
    of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
    The original matrix m stays intact. */
_const_  static
double *transpose(const double *, const unsigned , const unsigned , double *);

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant doesn't transpose matrix b, and it's a lot slower. */
_const_ static
double *dot_simple(const double *, const unsigned , const unsigned ,
                   const double *, const unsigned , const unsigned );

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant transposes matrix b, and it's a lot faster. */
_const_   static
double *dot(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, \
            const double *b, const unsigned n_rows_b, const unsigned n_cols_b);

/* Prints vector, or matrix. */
_const_ static
void print(const double *a, const unsigned n_rows_a, const unsigned n_cols_a);
/////////////////////////////////////////////////////////////////////////

/* Initializes vector or matrix, sequentially, with indices. */
_const_ static 
void init_seq(double *a, const unsigned n_rows_a, const unsigned n_cols_a) {
    for (unsigned i = 0; i < n_rows_a; i++) {
        for (unsigned j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = i*n_cols_a + j;
        }
    }
}

/* Initializes vector or matrix, randomly. */
_const_ static 
void init_rand(double *a, const unsigned n_rows_a, const unsigned n_cols_a) {
    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t j = 0; j < n_cols_a; j++) {
            a[i*n_cols_a + j] = rand() / (double)RAND_MAX;
        }
    }
}

/*  Takes and returns a new matrix, t, which is a transpose of the original one, m.
    It's also flat in memory, i.e., 1-D, but it should be looked at as a transpose
    of m, meaning, n_rows_t == n_cols_m, and n_cols_t == n_rows_m.
    The original matrix m stays intact. */
_const_ static 
double *transpose(const double *m, const unsigned n_rows_m, const unsigned n_cols_m, double *t) {
    for (size_t i = 0; i < n_rows_m; i++) {
        for (size_t j = 0; j < n_cols_m; j++) {
            t[j*n_rows_m + i] = m[i*n_cols_m + j];
        }
    }

    return t;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant doesn't transpose matrix b, and it's a lot slower. */
_const_ static
double *dot_simple(const double *a, const unsigned n_rows_a, const unsigned n_cols_a,
                   const double *b, const unsigned n_rows_b, const unsigned n_cols_b) {

    assert (n_cols_a == n_rows_b) ;

    double *c = 0 ; // malloc(n_rows_a * n_cols_b * sizeof(*c));

    ALLOC_WITH_POLICY( c, n_rows_a * n_cols_b * sizeof(*c));

    for (size_t i = 0; i < n_rows_a; i++) {
        for (size_t k = 0; k < n_cols_b; k++) {
            double sum = 0.0;
            for (size_t j = 0; j < n_cols_a; j++) {
                sum += a[i*n_cols_a + j] * b[j*n_cols_b + k];
            }
            c[i*n_cols_b + k] = sum;
        }
    }

    return c;
}

/* Dot product of two arrays, or matrix product
 * Allocates and returns an array.
 * This variant transposes matrix b, and it's a lot faster. */
_const_  static 
double *dot(const double *a, const unsigned n_rows_a, const unsigned n_cols_a, 
            const double *b, const unsigned n_rows_b, const unsigned n_cols_b) {

    assert (n_cols_a == n_rows_b) ;

    double * bt = 0 ;
    ALLOC_WITH_POLICY( bt, n_rows_b * n_cols_b * sizeof(*b) );

    double *c = 0;    
    ALLOC_WITH_POLICY( c, n_rows_a * n_cols_b * sizeof(*c) );


    bt = transpose(b, n_rows_b, n_cols_b, bt);

    for (unsigned i = 0; i < n_rows_a; i++) {
        for (unsigned k = 0; k < n_cols_b; k++) {
            double sum = 0.0;
            for (unsigned j = 0; j < n_cols_a; j++) {
                sum += a[i*n_cols_a + j] * bt[k*n_rows_b + j];
            }
            c[i*n_cols_b + k] = sum;
        }
    }
    return c;
}

#endif // 0

/* Prints vector, or matrix. */
_const_ static
void print(const double *a, const unsigned n_rows_a, const unsigned n_cols_a) {
    for (unsigned i = 0; i < n_rows_a; i++) {
        for (unsigned j = 0; j < n_cols_a; j++) {
            MATMUL_PRINT( "%8.3f ", a[i*n_cols_a + j]);
        }
        MATMUL_PRINT("\n");
    }
    MATMUL_PRINT("\n");
}


int main(int argc, char *argv[]) {
    _CRT_UNUSED(argc);
    _CRT_UNUSED( argv );
    /* Intializes random number generator */
    time_t t;
    srand((unsigned)time(&t));
    srand(0);

    /* For measuring time */
    double t0, t1;

    // const unsigned scale = 10; // provokes matrix print 
    const unsigned scale = 400;
    const unsigned n_rows_a = 4 * scale;
    const unsigned n_cols_a = 3 * scale;
    const unsigned n_rows_b = 3 * scale;
    const unsigned n_cols_b = 2 * scale;

    double * _cleanup_(cleanup_free)  a = NULL;
    double * _cleanup_(cleanup_free)  b = NULL;
    double * _cleanup_(cleanup_free)  c = NULL;
    double * _cleanup_(cleanup_free)  d = NULL;

    ALLOC_WITH_POLICY(a , n_rows_a * n_cols_a * sizeof(*a) );
    ALLOC_WITH_POLICY(b , n_rows_b * n_cols_b * sizeof(*b) );

    init_rand(a, n_rows_a, n_cols_a);
    init_rand(b, n_rows_b, n_cols_b);

    init_seq(a, n_rows_a, n_cols_a);
    init_seq(b, n_rows_b, n_cols_b);

    t0 = omp_get_wtime();
    c = dot_simple(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b);
    t1 = omp_get_wtime();
    MATMUL_PRINT("Dot Simple: Elapsed time %.3f s\n", t1 - t0);

    t0 = omp_get_wtime();
    d = dot(a, n_rows_a, n_cols_a, b, n_rows_b, n_cols_b);
    t1 = omp_get_wtime();
    MATMUL_PRINT("Dot: Elapsed time %.3f s\n", t1 - t0);

    if (scale < 11) {
        MATMUL_PRINT("Matrix A:\n");
        print(a, n_rows_a, n_cols_a);
        MATMUL_PRINT("Matrix B:\n");
        print(b, n_rows_b, n_cols_b);
        MATMUL_PRINT("Matrix C:\n");
        print(c, n_rows_a, n_cols_b);
        MATMUL_PRINT("Matrix D:\n");
        print(d, n_rows_a, n_cols_b);
    }

#ifdef WIN32
    system("pause");
#endif
    return EXIT_SUCCESS ;
}

