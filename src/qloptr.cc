/**
 * @file qlopt.cc
 *
 * @date 24.06.2016
 * @author: M. Baaske
 */

#include "qsoptim.h"
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

/* R Interface functions  */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, (n)}

R_NativePrimitiveArgType myC_t[] = {
	REALSXP, INTSXP, INTSXP
};

static R_CMethodDef CEntries[] = {
	{"isPositiveDefinite", (DL_FUNC) &isPositiveDefinite, 3, myC_t},
	{ NULL, NULL, 0, NULL}
};

static R_FortranMethodDef FEntries[] = {
	{"dggev", (DL_FUNC) &F77_NAME(dggev), 17, NULL},
	{ NULL, NULL, 0, NULL}
};

static R_CallMethodDef CallEntries[] = {
    CALLDEF(kriging,5),
    CALLDEF(estimateJacobian,5),
    CALLDEF(getDualKrigingWeights,3),
    CALLDEF(Fmat,2),
    CALLDEF(Pmat,1),
    CALLDEF(QSopt,7),
	CALLDEF(initQL,5),
    CALLDEF(finalizeQL,0),
    CALLDEF(qDValue,1),
	CALLDEF(mahalValue,1),
	CALLDEF(mahalanobis,7),
	CALLDEF(quasiDeviance,7),
    CALLDEF(covMatrix,2),
	CALLDEF(covValue,2),
    {NULL, NULL, 0}
};

#ifdef  __cplusplus
extern "C" {
#endif

void R_init_qle(DllInfo *info)
{
  R_registerRoutines(info, CEntries, CallEntries, FEntries, NULL);
  R_useDynamicSymbols(info, FALSE);
}


/*
 void R_unload_qle(DllInfo *info){
  /// Release resources.
}
*/


#ifdef  __cplusplus
}
#endif



