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

#define USE_R
#ifndef DEBUG_H_
#define DEBUG_H_

#include "basics.h"

#define NO_ERROR 0
#define NO_WARNING 0

#define MAX_ERRORS 100
#define MAX_WARNINGS 100
#define maxErrorStr 1000
#define nErrorLocations 100
#define nWarningLocations 100

#define ERROR_NUM_MIN 100
#define ERROR_NUM_MAX 500
#define WARNING_NUM_MIN 101

#define BEGIN_ERR { Rprintf("\n\n =====================  ERROR block ========================\n"); }
#define BEGIN_WRR { Rprintf("\n\n =====================  WARNING block ======================\n"); }

#define PRINT_MSG(s)  Rprintf("[%s:%u]%s\n", __FILE__, (unsigned)__LINE__, s)

extern int  PL;
extern char C_MSG_BUFFER[1000];
extern char R_MSG_BUFFER[1000];

extern char ERROR_STRING[maxErrorStr];
extern char ERROR_LOC[nErrorLocations];


#define LOG_ERROR(X, INFO, MSG) { \
  errorMSG(X, C_MSG_BUFFER); \
  std::sprintf(R_MSG_BUFFER, "Error message: `%s', error code %d: %s",C_MSG_BUFFER, X,MSG); \
  setError(X,R_MSG_BUFFER,__LINE__, __FILE__, INFO); \
}

#define LOG_WARNING(X, INFO, MSG) { \
  warningMSG(X, C_MSG_BUFFER); \
  std::sprintf(R_MSG_BUFFER, "Warning message: `%s', warning code %d: %s",C_MSG_BUFFER, X, MSG); \
  setWarning(X,R_MSG_BUFFER,__LINE__, __FILE__, INFO); \
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


#ifdef USE_R
#else
  #include <iostream>
  #include <sstream>
  #include <stdio.h>
#endif


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

typedef enum  warning_type_s {
      GENERIC_WARNING = 1,
      FINAL_WARNING = 1000
} warning_type;


typedef struct errorStatus_s {
  char msg[maxErrorStr];       /* some additional message */
  char file[nErrorLocations]; /* the file */

  int id,       /* error id see error_type */
      code,     /* specific error code  */
      numErr,   /* current error number */
      line,     /* the line */
      isError;  /* error = 1, warning = 0 */

} errorStatus;

typedef struct warningStatus_s {
  char msg[maxErrorStr];       /* some additional message */
  char file[nErrorLocations]; /* the file */

   int id,              /* error id see error_type */
      code,             /* specific error code  */
     numWrr,            /* current error number */
       line,            /* the line */
       isWarning;       /* error = 1, warning = 0 */

} warningStatus;

typedef warningStatus *warning_queue[MAX_WARNINGS];
typedef errorStatus *error_queue[MAX_ERRORS];

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
extern void setError(int, const char *, int ,const char *, int);
extern void setWarning(int , const char* , int , const char* , int);
extern void printMSG(int, const char *, int, const char * );

extern void cleanErrors();
extern void incWarning();
extern void incError();
extern void cleanWarnings();

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

#endif /* DEBUG_H */
