/**
 * @file        qloptr.h
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface)
 *
 * Explanation:
 *
 */

#ifndef SYSKRIGE_H_
#define SYSKRIGE_H_

#include "basics.h"
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

#ifdef  __cplusplus
extern "C" {
#endif

SEXP emptyErrorCache();
SEXP emptyWarningCache();

/* DLL Init */
void R_init_qle(DllInfo *info);

#ifdef  __cplusplus
}
#endif

#endif /* SYSKRIGE_H_ */

