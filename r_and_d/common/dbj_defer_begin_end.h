#ifndef DBJ_DEFER_BEGIN_END_H
#define DBJ_DEFER_BEGIN_END_H
/*
   https://godbolt.org/z/zsYre5Kfj

  first seen it here https://youtu.be/QpAhX-gsHMs
  2021-APR  dbj@dbj.org added macro defer()

  Two Macros. Usage:
  
  // these two footprints are totaly whatever you need!
  // Q: when are the macro arguments executed? 
  // A: be sure to understand ;)
  void on_scope_begin(void) ;
  void on_scope_end(void) ;

  beginend( on_scope_begin(), on_scope_end() ) { }

  defer( on_scope_end() ) { }

  See the ad-hoc-demo bellow.

  Q: What might be also good about this two macros? 
  A: They work on any C version, using any compiler.

  There is no trick: Both macros are for(){ } that loops once
  I see no reason for this not to work on any C
*/
#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>

#define macro_concat_(a,b) a ## b
#define macro_concat(a,b) macro_concat_(a,b)
#define macro_var(name) macro_concat( name, __LINE__)

#define beginend( start,end) \
for ( int macro_var(_i_) = ( start, 0 ) ; ! macro_var(_i_) ; \
( macro_var(_i_) += 1, end) )

static inline void do_nothing_ (void) {}
// #define defer( at_scope_end ) beginend( __LINE__, at_scope_end ) 
#define defer( at_scope_end) for ( int macro_var(_i_) = 0;\
 ! macro_var(_i_) ; ( macro_var(_i_) += 1, at_scope_end) )

/*
 Ad-hoc demo
*/
#if DBJ_DEFER_BEGIN_END_DEMO

#include <string.h>
#include <stdio.h>
#include <stdbool.h>

// NOTE: F must be a string literal
#define B(X_) (X_) ? "true" : "false"
#define P(F,X_) fprintf( stdout, "\n%04d : %16s :\t" F, __LINE__, #X_, (X_))
#define M(S) fprintf( stdout, "\n%04d : %16s :\t%s", __LINE__, " ", S)

#ifdef __linux__
#include <unistd.h> // readlink
#endif

static void here_ () {  M("Begin!");  }
static void there_ () {  M("End!");  }

static void make_ ( char ** ptr, size_t sze_ ) {  *ptr= calloc(sze_, sizeof(char)); printf("\nmade  (%p) char [%zu]", *ptr , sze_); }
static void free_ ( void * ptr ) {printf("\nfreed (%p)", ptr ); free(ptr); ptr = NULL;  }

static void close_file ( FILE * fp_ ) {
    if ( fp_ == NULL) {  M("close_file(): Avoiding NULL file pointer"); 
    return ; }
    if ( ferror(fp_))  perror("I/O error");
    P("closing FILE * (%p)", (void*)fp_);
    fclose(fp_ ) ;  
 }

int main (void)
{
    beginend( here_(), there_() )
    {
         M("Work.");
    }

    char * bufy_ = 0 ;
    printf("\n");

    beginend( make_( &bufy_, 0xFF ), free_( bufy_ ) )
    {
         memcpy( bufy_, "DATA", 5 );
         printf("\nusing (%p) : \"%s\"", (void*)bufy_, bufy_);
    }

    printf("\n");
    FILE * dummsy_ = tmpfile();

     defer( close_file( dummsy_ )) {
         P("temp FILE * (%p)", (void*)dummsy_);

#ifdef __linux__
    // Linux-specific 
    char fname[FILENAME_MAX] = {0}, link[FILENAME_MAX] = {0};
    sprintf(fname, "/proc/self/fd/%d", fileno(dummsy_));
        P("temp FILE name: %s", fname);
    if(readlink(fname, link, sizeof link - 1) > 0)
        P("temp FILE link: %s", link);
#endif // __linux__

    } 

    return 42 ;
}

#undef  B
#undef  P
#undef  M

#endif // DBJ_DEFER_BEGIN_END_DEMO

#undef macro_concat_
#undef macro_concat
#undef macro_var

#endif // DBJ_DEFER_BEGIN_END_H