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

#include <sstream>
#include <cstring>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("qle", String)
#else
#define _(String) (String)
#endif

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

#endif



