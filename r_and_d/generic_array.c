/*
(c) 2021 by dbj@dbj.org -- LICENSE_DBJ -- https://dbj.org/license_dbj/
*/
#include <stdlib.h>
#include <assert.h>
#include <string.h>  // memcpy

#define _cleanup_(x) __attribute__((__cleanup__(x)))

// usage: type_name _auto_cleanup_ pointer_name ;
// example
// char * _auto_cleanup_ str = strdup("string content");
//
#define _auto_cleanup_ _cleanup_(cleanup_free)

// DBJ 2021-JUL-05 added
static inline void cleanup_free(void* p) { free(*(void**)p); }

#define CCAT_(a, b) a##b
#define CCAT(a, b) CCAT_(a, b)

// be aware you can assign struct instance to instance
// but if flexy array is inside it will not be copied over

#define ARRTYPE(T_)                                                          \
  typedef struct CCAT(T_, S_) {                                              \
    size_t size;                                                             \
    T_ data[];                                                               \
  } CCAT(T_, _array);                                                        \
                                                                             \
  static inline CCAT(T_, _array) *                                           \
      CCAT(new_, CCAT(T_, _array))(unsigned size_) {                         \
    CCAT(T_, _array)* retval = calloc(1, sizeof(CCAT(T_, _array)) + size_);  \
    if (retval) retval->size = size_;                                        \
    return retval;                                                           \
  }                                                                          \
                                                                             \
  static inline CCAT(T_, _array) *                                           \
      CCAT(copy_, CCAT(T_, _array))(CCAT(T_, _array) * src_) {               \
    assert(src_);                                                            \
    CCAT(T_, _array)* retval =                                               \
        calloc(1, sizeof(CCAT(T_, _array)) + src_->size);                    \
    if (retval) memcpy(retval, src_, sizeof(CCAT(T_, _array)) + src_->size); \
    return retval;                                                           \
  }

// declares and defines 
// float_array
// new_float_array( unsigned );
// copy_float_array( unsigned );
ARRTYPE(float);

// char_array
// new_char_array( unsigned );
// copy_char_array( unsigned );
ARRTYPE(char);

int main(void) 
{
  float_array* _auto_cleanup_ f255 = new_float_array(0xFF);
  f255->data[0] = 3.14;
  float_array* _auto_cleanup_ f255copy = copy_float_array(f255);
  assert( f255->data[0] == f255copy->data[0]);
  //
  char_array* _auto_cleanup_ c255 = new_char_array(0xFF);
  c255->data[0] = '!';
  char_array* _auto_cleanup_ c255copy = copy_char_array(c255);
  assert( c255->data[0] == c255copy->data[0] );

  return 42;
}

  