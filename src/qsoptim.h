/**
 * @file        qsoptim.h
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions for optimization
 *
 * Explanation: Fisher Quasi-Scoring iteration based on surrogate statistics
 */

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

#include "kriging.h"

typedef enum {
	 QFS_ERROR = -10, 				/* generic failure code */
	 QFS_MAXITER_REACHED = -5,
     QFS_LINESEARCH_FAILURE = -4,
     QFS_LINESEARCH_ZEROSLOPE = -3,
     QFS_BAD_DIRECTION = -2,    	/* Calculating Newton direction failed*/
	 QFS_NO_CONVERGENCE = -1,   	/* no convergence */
     QFS_CONVERGENCE = 1,       	/* convergence */
     QFS_SCORETOL_REACHED = 2,
     QFS_FTOLREL_REACHED = 3,
     QFS_STOPVAL_REACHED = 4,
     QFS_XTOL_REACHED = 5,
     QFS_GRADTOL_REACHED = 6,
     QFS_SLOPETOL_REACHED = 7,
     QFS_LOCAL_CONVERGENCE = 10
} qfs_result;

typedef struct qfs_options_s {
  ql_model qlm;

  int num_iter, num_eval; /* used */
  int pl, info;	 		  /* print level */

  double grad_tol,        /* stopping criteria */
  	     ftol_stop,
		 ftol_abs,
		 ftol_rel,
		 score_tol,
		 xtol_rel,
		 slope_tol;

  int max_iter;    /* limits */

  qfs_options_s(ql_model _qlm, SEXP R_options) : qlm(_qlm), num_iter(0), num_eval(0)
  {
	info = 0;
    pl = asInteger(getListElement( R_options, "pl"));

	ftol_rel  = asReal(getListElement( R_options, "ftol_rel" ));
	ftol_stop = asReal(getListElement( R_options, "ftol_stop"));
	ftol_abs  = asReal(getListElement( R_options, "ftol_abs"));
	score_tol = asReal(getListElement( R_options, "score_tol"));
	xtol_rel  = asReal(getListElement( R_options, "xtol_rel" ));
	grad_tol  = asReal(getListElement( R_options, "grad_tol" ));
	slope_tol = asReal(getListElement( R_options, "slope_tol"));
	max_iter  = asInteger(getListElement( R_options, "maxiter"));
  }

} qfs_options_t, *qfs_options;

#ifdef  __cplusplus
extern "C" {
#endif

SEXP emptyErrorCache();
SEXP emptyWarningCache();

// Quasi-Scoring iteration
SEXP QSopt(SEXP R_start, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm, SEXP R_opt);

#ifdef  __cplusplus
}
#endif

#endif /* AUX_H_ */

