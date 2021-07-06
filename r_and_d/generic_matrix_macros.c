#if ! defined(__clang__) 
// #error clang required
#endif // ! __clang__

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

/*
this is developed in VS 2019 using clang 11.0
all is fine beside intellisence going berserk ... sometimes
*/


#define CONCAT_IMPL( x, y ) x##y
#define MACRO_CONCAT( x, y ) CONCAT_IMPL( x, y )

/*
 DBJ_MX_STRUCT(int,3,4) int_matrix ;  makes

 struct { 
    const size_t width; 
    const size_t height; 
    int(*data)[3][4];
} mx ;

  * mx.data        -> type is int[w][h]
    mx.data[0][0]  -> type is int[w]
  * mx.data[0][0]  -> type is int
& * mx.data[0][0]  -> type is int *

NOTE! above mx fields are totally uninitialized, be carefull.

NOTE! if const is used bellow one can not swap two matrixes with DBJ_SWAP()

*/
#define DBJ_MX_STRUCT(T,WIDTH,HEIGHT) \
struct { \
    /*const*/ size_t width; \
    /*const*/ size_t height; \
    T(*data)[WIDTH][HEIGHT]; \
} 
/*
dbj note: above could be T data[WIDTH][HEIGHT] ... 
which would likely deplete your stack rather quickly
but yes, this is how embeded OS does it ... all on a  stack, no heap
*/

/*
matrix element pointer 
move these asserts inside a macro if required:
    assert( MX_.data ); 
    assert( R_ < MX_.width ); 
    assert( C_ < MX_.height); 
*/
#define DBJ_MX_EP(MX_, R_, C_) (& (*MX_.data)[R_][C_]) 

/*
https://www.ibm.com/support/knowledgecenter/en/ssw_ibm_i_71/rzarg/typeof_operator.htm

dbj: bellow is totaly generic. if used for matrix as declared above
     one could check that sizes are equal
     bellow will fail if a and b types are sufficiently different
     but if they are not, swap will happen and heisenbug is born
*/
#define DBJ_SWAP(a,b) do { __typeof__(a) temp; temp = a; a = b; b = temp; } while(0)

#define DBJ_MX_TYPE(mx_) __typeof__(mx_) 
// remember this yields type T *, not  T 
#define DBJ_MX_EP_TYPE(mx_) __typeof__( & (*mx_.data)[0][0] ) 
// this is element type, aka T
#define DBJ_MX_E_TYPE(mx_) __typeof__( (*mx_.data)[0][0] ) 

// make mx instance as compund literal
// us this ad-hoc and localy 
// NOTE: calloc/malloc can not be used from a global space
#define DBJ_MX_INIT(T, MXT_,W_,H_) \
(MXT_){ W_, H_, calloc(1, sizeof(T[W_][H_])) }

// make a mx typedef
// also make a factory function creating the instance of
// factory. Function name is TYPENAME_##_factory
// example: 
// my_matrix_type  my_matrix_type_factory() ;
// returns instances of my_matrix_type

#define DBJ_MX_TYPE_AND_FACTORY(T_, W_, H_, TYPENAME_ ) \
    typedef DBJ_MX_STRUCT(T_, W_, H_) TYPENAME_ ; \
\
    static inline TYPENAME_ TYPENAME_##_factory () { \
        return DBJ_MX_INIT(T_, TYPENAME_,W_,H_) ; \
    }

/*
bellow is an callback as a macro 
(example) : zero the element of the matrix
If flag is set to false, for each stops

move these two asserts inside a macro if required:
    assert( ROW_ < MX_.width ); \
    assert( COL_ < MX_.height); \
*/
#define callback_zoro(ROW_,COL_, MX_, BOOL_PTR_FLAG_ ) \
do { \
    (*MX_.data)[ROW_][COL_] = ( DBJ_MX_E_TYPE(MX_) )0; \
    *BOOL_PTR_FLAG_ = true; \
} while(0)

/*
try to use it all through for-each
note: CBACK can be macro too
*/
#define DBJ_MX_FOREACH_( MX_,  CBACK, GOTOLABEL ) \
do { \
    for (size_t i = 0; i < MX_.width ; i++) \
        for (size_t j = 0; j < MX_.height; j++) { \
            bool retval_ = true; \
            CBACK(i, j, MX_, &retval_ ); \
            if (!retval_) goto GOTOLABEL ;\
        } \
    GOTOLABEL : ;\
} while(0)


/*
here we create a unique goto label
so we can repeat this macro inside a same scope
*/
#define DBJ_MX_FOREACH( MX_,  CBACK ) \
DBJ_MX_FOREACH_( MX_,  CBACK, MACRO_CONCAT( dbj_exit_label_, __COUNTER__) ) 

#define DBJ_MX_PRINT(FMT_, MX_ )\
do {\
    for ( size_t i = 0; i < MX_.width ; i++) {\
        for (size_t j = 0; j < MX_.height ; j++) {\
            printf(FMT_, *MX_.data[i][j]);\
        }\
        printf("\n"); \
    }\
} while(0)

//
// note: this is naive algorithm
//
// standard declaration
//static void multiply(int m, int n, int p, int a[m][n], int b[n][p], int c[m][p]);
//
// using DBJ MX types
//static void multiply(int_3_4_mx a, int_4_3_mx b, int_3_3_mx c )
//{
//    const size_t m = a.width, n = a.height, p = b.width;
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < p; j++) {
//            * c.data[i][j] = 0;
//            for (int k = 0; k < n; k++) {
//                (* c.data[i][j]) += (* a.data[i][k]) * (* b.data[k][j]);
//            }
//        }
//    }
//}

/*
core macros test
to see this working (or not) follow carefully through debugger
*/
static void dbj_mx_core_tests()
{
    enum { W = 2, H = 2 };
    // define and init the matrix variable
    // NOTE: using VLA in a sizeof 
    // NOTE: can not call calloc/malloc from a global space!
    DBJ_MX_STRUCT(int, W, H) mx = { W, H, calloc(1, sizeof(int[W][H])) };

    DBJ_MX_TYPE(mx) mx2 = { W, H, calloc(1, sizeof(int[W][H])) };

    DBJ_MX_E_TYPE(mx)   integer = 42;
    DBJ_MX_EP_TYPE(mx)  integer_ptr = (DBJ_MX_EP_TYPE(mx)) & integer;

    /* due to usage of calloc all is zeroed , but we want to test
    callback_zoro to  set the element 1,1 to 0 */
    bool callback_flag;
    callback_zoro(1, 1, mx, &callback_flag);

    /* give some values */

    *DBJ_MX_EP(mx, 0, 0) = 1;
    *DBJ_MX_EP(mx, 0, 1) = 2;
    *DBJ_MX_EP(mx, 1, 0) = 3;
    *DBJ_MX_EP(mx, 1, 1) = 4;

    *DBJ_MX_EP(mx2, 0, 0) = 5;
    *DBJ_MX_EP(mx2, 0, 1) = 6;
    *DBJ_MX_EP(mx2, 1, 0) = 7;
    *DBJ_MX_EP(mx2, 1, 1) = 8;

    /* we can swap two matrixes of equal type */
    DBJ_SWAP(mx, mx2);
}

/*
Another test bench
*/
enum { r1 = 3, c1 = 4, r2 = c1, c2 = r1 };

// make matrix types declarations 
// and factory functions implementations
DBJ_MX_TYPE_AND_FACTORY(int,r1,c1, int34mx );
DBJ_MX_TYPE_AND_FACTORY(int,r2,c2, int43mx );
DBJ_MX_TYPE_AND_FACTORY(int,r1,c2, int33mx );

int dbj_matrix_struct_test ()
{
    dbj_mx_core_tests();

    int34mx mx34 = int34mx_factory();
    int43mx mx43 = int43mx_factory();
    int33mx mx33 = int33mx_factory();

    // callback_zoro is also a macro
    DBJ_MX_FOREACH( mx34, callback_zoro );
    DBJ_MX_FOREACH( mx43, callback_zoro );
    DBJ_MX_FOREACH( mx33, callback_zoro );

    printf("\nmatrix mx34\n");
    DBJ_MX_PRINT("%4d ", mx34);
    printf("\nmatrix mx43\n");
    DBJ_MX_PRINT("%4d ", mx43);
    printf("\nmatrix mx33\n");
    DBJ_MX_PRINT("%4d ", mx33);

    return 42;
}

int main(void)
{
    return dbj_matrix_struct_test() ;
}

/*
undefine macros that can clash with other macros
*/
#undef CONCAT_IMPL
#undef MACRO_CONCAT

