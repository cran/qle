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

#define maxErrorStr 100
#define nErrorLocations 100

#define ERROR_NUM_MIN 100
#define ERROR_NUM_MAX 500
#define WARNING_NUM_MIN 101

#define BEGIN_ERR { Rprintf("\n\n =====================  ERROR block ========================\n"); }
#define BEGIN_WRR { Rprintf("\n\n =====================  WARNING block ======================\n"); }

#define PRINT_MSG(s)  Rprintf("[%s:%u]%s\n", __FILE__, (unsigned)__LINE__, s)

#define PRINT Rprintf
#define PRINTF Rprintf

extern int  PL;
extern char C_MSG_BUFFER[100];
extern char R_MSG_BUFFER[100];
extern char ERROR_LOC[nErrorLocations];


#define LOG_ERROR(X, INFO, MSG) { \
  errorMSG(X, C_MSG_BUFFER); \
  std::sprintf(R_MSG_BUFFER, "Error: `%s', error code %d: %s", C_MSG_BUFFER, X, MSG); \
  PRINT(" %s ", R_MSG_BUFFER); \
}

#define LOG_WARNING(X, INFO, MSG) { \
  warningMSG(X, C_MSG_BUFFER); \
  std::sprintf(R_MSG_BUFFER, "Warning: `%s', warning code %d: %s", C_MSG_BUFFER, X, MSG); \
  PRINT(" %s ", R_MSG_BUFFER); \
}


/* trigger error and warnings */
#define ERR(X) { \
    std::sprintf(ERROR_LOC, "in `%s', at %u\n", __FILE__, (unsigned)__LINE__); \
    std::sprintf(R_MSG_BUFFER, "%s, %s\n",ERROR_LOC, X); \
    error(_(R_MSG_BUFFER)); \
}

#define XERR(X,MSG) { \
     BEGIN_ERR; \
     errorMSG(X, C_MSG_BUFFER); \
     std::sprintf(ERROR_LOC, "in `%s', at %u\n", __FILE__, (unsigned)__LINE__); \
     std::sprintf(R_MSG_BUFFER, "%s, error code %d. (%s) %s", ERROR_LOC, X, MSG, C_MSG_BUFFER); \
     error(_(R_MSG_BUFFER)); \
  }

// X is message string
#define WRR(X) { \
    std::sprintf(ERROR_LOC, "in `%s', at %u\n", __FILE__, (unsigned)__LINE__); \
    std::sprintf(R_MSG_BUFFER, "%s, %s\n",ERROR_LOC, X); \
    warning(_(R_MSG_BUFFER)); \
}

// X is integer error code
#define XWRR(X,MSG) { \
    BEGIN_WRR; \
    errorMSG(X, C_MSG_BUFFER); \
    std::sprintf(ERROR_LOC, "in `%s', at %u\n", __FILE__, (unsigned)__LINE__); \
    std::sprintf(R_MSG_BUFFER, "%s, warning code %d: (%s) %s", ERROR_LOC, X, MSG, C_MSG_BUFFER); \
    warning(_(R_MSG_BUFFER)); \
}


#define PERR(X) { \
     BEGIN_ERR; \
     std::sprintf(R_MSG_BUFFER, "%s\n%s: %s", ERROR_LOC, param, X); error(R_MSG_BUFFER); \
}

// memory allocation errors
#define MEM_ERR(n, t) { \
  std::sprintf(ERROR_LOC, "'calloc' error in `%s', at %u\n", __FILE__, (unsigned)__LINE__); \
  std::sprintf(R_MSG_BUFFER, "%s, %s (%.0f of %u bytes) \n",ERROR_LOC, "Could not allocate memory.", (double) (n), (unsigned)sizeof(t)); \
  error(_(R_MSG_BUFFER)); \
}

typedef enum  {
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
        PARAM_ERROR = 500,
        FINAL_ERRROR = 1000      /*not to be changed */

} error_type;

template<class Type>
void printArray(const char fmt[], const Type *v, int *lx) {
  Rprintf("[");
  for(int i=0; i < *lx; i++)
     Rprintf( fmt , v[i]);
  Rprintf("]");
}


/////////////////////// Error /////////////////////////////////////////////////////////////////////////////////////
extern void errorMSG(int, char*);
extern void warningMSG(int , char *, char* );
extern void printMSG(int, const char *, int, const char * );

void printMatrix(const char ch[], const double *mat, int *irows, int *icols);
void printVector(const char ch[], const double *vec, int *lx);

void print_R_matrix( SEXP A, const std::string& matrix_name);
void print_R_vector( SEXP v, const std::string& vector_name);

#define DEBUG_INFO(_x) \
	do { \
		fprintf(stderr, "===== BEGIN: DEBUG block =====\n" \
				"file: %s, line: %u\n", \
				__FILE__, (unsigned)__LINE__); \
				_x; \
				fprintf(stderr, "===== E N D: DEBUG block =====\n"); \
				fflush(stderr); \
	} while (0)

#define DEBUG_PRINT_VECTOR_R( c , cstr ) { print_R_vector( c , cstr ); }
#define DEBUG_DUMP_VAR(x,fmt) { PRINT("%s:%u: %s=\n" fmt, __FILE__, (unsigned)__LINE__, #x, x); }
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
