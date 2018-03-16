/**
 * @file        error.h
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       Debugging definitions and error functions
 *
 *
 * Explanation: Trigger errors and warnings either directly by calling
 *              ERR, WRR, XERR, XWRR makros or store them in some
 *              objects to read from R and trigger them as you want
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#include "basics.h"

#include <R.h>
#include <Rdefines.h>

#define NO_ERROR 0
#define NO_WARNING 0
#define INTERNAL_ERROR 10

extern int  PL;
extern char C_MSG_BUFFER[1000];

#define PRINT_MSG(M) { \
  Rprintf(" %s \n %s (line=%u)\n", M,  __FILE__, (unsigned)__LINE__); \
}

#define LOG_ERROR(X, NAME) { \
  errorMSG(X, NAME); \
  if(PL > 0 ) \
  	Rprintf("%s ", C_MSG_BUFFER); \
}

#define LOG_WARNING(X, NAME) { \
  warningMSG(X, NAME); \
  if(PL > 0 ) \
   	Rprintf("%s ", C_MSG_BUFFER); \
}


/* trigger error and warnings */
#define ERR(MSG) { \
  Rprintf("error at %s (line=%u)\n", __FILE__, (unsigned)__LINE__); \
  Rf_error(_(MSG)); \
}

#define XERR(X,NAME) { \
     errorMSG(X, NAME); \
     Rprintf("%s  (line=%u)\n", __FILE__, (unsigned)__LINE__); \
     Rf_error(_(C_MSG_BUFFER)); \
  }

// X is message string
#define WRR(MSG) { \
  Rprintf("warning at %s (line=%u)\n", __FILE__, (unsigned)__LINE__); \
  Rf_warning(_(MSG)); \
}

// X is integer error code
#define XWRR(X,NAME) { \
    warningMSG(X, NAME); \
    Rprintf("%s (line=%u)\n", __FILE__, (unsigned)__LINE__); \
    Rf_warning(_(C_MSG_BUFFER)); \
}

// memory allocation errors
#define MEM_ERR(n, t) { \
  Rprintf("%s (line=%u)\n", __FILE__, (unsigned)__LINE__); \
  Rprintf("(%.0f of %u bytes) \n", (double) (n), (unsigned)sizeof(t)); \
  Rf_error(_("Could not allocate memory.")); \
}

typedef enum {
        SYSTEM_ERROR = 1,
        ALLOC_ERROR = 2,
        MEMORY_ERROR = 3,
        MATH_ERROR = 4,
        NaN_ERROR = 5,
	    LAPACK_ERROR = 100,
        LAPACK_QR_ERROR = 101,
        LAPACK_PMAT_ERROR = 102,
        LAPACK_SOLVE_ERROR = 103,
        LAPACK_INVERSION_ERROR = 104,
        LAPACK_FACTORIZE_ERROR = 105,
        FINAL_ERRROR = 1000      /*not to be changed */

} error_type;

typedef enum {
        NaN_WARNING = 1,
		POSDEF_WARNING = 2,
		LAPACK_WARNING = 3,
		GENERIC_WARNING = 500,
		FINAL_WARNING = 1000      /*not to be changed */
} warning_type;

template<class Type>
void printArray(const char fmt[], const Type *v, int *lx) {
  Rprintf("[");
  for(int i=0; i < *lx; i++)
     Rprintf( fmt , v[i]);
  Rprintf("]");
}


/////////////////////// Error /////////////////////////////////////////////////////////////////////////////////////
extern void errorMSG(int, const char*);
extern void warningMSG(int , const char* );

void printMatrix(const char ch[], const double *mat, int *irows, int *icols);
void printVector(const char ch[], const double *vec, int *lx);

void print_R_matrix( SEXP A, const std::string& matrix_name);
void print_R_vector( SEXP v, const std::string& vector_name);

#define DEBUG_INFO(_x) \
	do { \
		fprintf(stderr, "===== B E G I N: DEBUG block =====\n" \
				"file: %s, line: %u\n", \
				__FILE__, (unsigned)__LINE__); \
				_x; \
				fprintf(stderr, "===== E N D: DEBUG block =====\n"); \
				fflush(stderr); \
	} while (0)

#define DEBUG_PRINT_VECTOR_R( c , cstr ) { print_R_vector( c , cstr ); }
#define DEBUG_DUMP_VAR(x,fmt) { Rprintf("%s:%u: %s=\n" fmt, __FILE__, (unsigned)__LINE__, #x, x); }
#define DEBUG_PRINT_MATRIX_R( C , cstr) { print_R_matrix( C , cstr); }

/////////////////// Some R convenience macros ////////////////////////////////////

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define GET_DIMS(A) INTEGER(coerceVector (getAttrib( (A), (R_DimSymbol) ) , INTSXP) )

#define SET_CLASS_NAME(RObject,ClassName) {            \
  SEXP RClass;                                         \
  PROTECT(RClass=allocVector(STRSXP,1));               \
  SET_STRING_ELT( (RClass) ,0, mkChar( (ClassName) ));  \
  classgets( (RObject), (RClass) );                    \
  UNPROTECT(1);                                        \
}

// Dimension names of a matrix, UNPROTECT later
#define SET_DIMNAMES_MATRIX(RObject,RNamedObject){ 		\
  SEXP R_names = getAttrib(RNamedObject,R_NamesSymbol); \
  SET_VECTOR_ELT(RObject, 0, R_names);					\
  SET_VECTOR_ELT(RObject, 1, R_names);					\
}

// Dimension names of a matrix, UNPROTECT later
#define SET_DIMNAMES_MATRIX2(RObject,RNamedObject0,RNamedObject1){ 		\
  SEXP R_names0 = getAttrib(RNamedObject0,R_NamesSymbol); 				\
  SEXP R_names1 = getAttrib(RNamedObject1,R_NamesSymbol); 				\
  SET_VECTOR_ELT(RObject, 0, R_names0);									\
  SET_VECTOR_ELT(RObject, 1, R_names1);									\
}


#endif /* DEBUG_H */
