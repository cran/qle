/**
 * @file        basics.h
 * @author      M. Baaske
 * @date        27/12/2012
 * @brief       Includes and basic macros
 *
 *
 */

#ifndef R_BASICS_H
#define R_BASICS_H 1

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <iostream>
#include <sstream>
#include <cstring>

#include <stdlib.h>
#include <assert.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("qle", String)
#else
#define _(String) (String)
#endif

/////////////////// Some R convenience macros ////////////////////////////////////

#define PRINT Rprintf
#define PRINTF Rprintf

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define GET_DIMS(A) INTEGER(coerceVector (getAttrib( (A), (R_DimSymbol) ) , INTSXP) )

// MEMORY ------------------------------------------------------------------------
#define FREE(p)  if ( (p) != NULL ) { std::free ( (void *) (p) ); (p) = NULL; }
#define DELETE(p) if ( (p) != NULL ) { delete (p); (p) = NULL;}
#define DELETE_ARR(p) if ( (p) != NULL ) { delete[] (p); (p) = NULL;}
#define MEMZERO(p, n) std::memset( (p), 0, (size_t)(n) * sizeof(*p) )
#define MEMCPY(p, q, n) std::memcpy( (p) , (q), (size_t)(n) * sizeof(*p) )
#define MALLOC(n, t) (t*) std::malloc( (size_t) (n), sizeof(t) )
#define CALLOC(n, t) (t*) std::calloc( (size_t) (n), sizeof(t) )


#define CALLOCX(p, n, t) { \
	if( ( (p) = CALLOC(n, t) ) == NULL ) \
	  MEM_ERR( (size_t) (n), t) \
}

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

//////////////////////////////////////////////////////////////////////////////////

#endif



