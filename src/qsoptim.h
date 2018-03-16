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
	 QFS_EVAL_ERROR = -8,			/* monitor function evaluation */
	 QFS_STEPTOL_REACHED = -6,		/* line search could not find sufficient decrease */
	 QFS_LINESEARCH_ERROR = -5,
	 QFS_MAXITER_REACHED = -4,
     QFS_LINESEARCH_ZEROSLOPE = -3,
     QFS_BAD_DIRECTION = -2,    	/* Calculating Newton direction failed*/
	 QFS_NO_CONVERGENCE = -1,   	/* no convergence */
     QFS_CONVERGENCE = 0,       	/* convergence */
     QFS_SCORETOL_REACHED = 1,
     QFS_FTOLREL_REACHED = 2,
     QFS_STOPVAL_REACHED = 3,
     QFS_SLOPETOL_REACHED = 4,
	 QFS_LOCAL_CONVERGENCE = 5,		/* approximate stationary point found at (scaled) norm^2 of quasi-score */
	 QFS_STEPMIN_REACHED = 6,		/* relative length of Newton direction is near zero, usually indicates convergence  */
	 QFS_XTOL_REACHED = 10			/* might have converged, however, sometimes also indicates problems at bounds */
} qfs_result;

typedef struct qfs_options_s {
  ql_model qlm;

  int num_iter, num_eval; /* used */
  int pl, info, doIobs;	  /* print level */

  double grad_tol,        /* stopping criteria */
  	     ftol_stop,
		 ftol_abs,
		 ftol_rel,
		 score_tol,
		 slope_tol,
		 xtol_rel;

  double *typf, *typx;

  int max_iter;    /* limits */

  qfs_options_s(ql_model _qlm, SEXP R_options) :
	  qlm(_qlm), num_iter(0), num_eval(0), pl(0), info(0), doIobs(FALSE),
	  typf(0), typx(0)
  {
    pl = asInteger(getListElement( R_options, "pl"));
    doIobs = asInteger(getListElement( R_options, "Iobs"));
	ftol_rel  = asReal(getListElement( R_options, "ftol_rel" ));
	ftol_stop = asReal(getListElement( R_options, "ftol_stop"));
	ftol_abs  = asReal(getListElement( R_options, "ftol_abs"));
	score_tol = asReal(getListElement( R_options, "score_tol"));
	xtol_rel  = asReal(getListElement( R_options, "xtol_rel" ));
	slope_tol = asReal(getListElement( R_options, "slope_tol" ));
	grad_tol  = asReal(getListElement( R_options, "grad_tol" ));
	max_iter  = asInteger(getListElement( R_options, "maxiter"));

	// scaling vectors
	typf = REAL(getListElement( R_options, "fscale" ));
	typx = REAL(getListElement( R_options, "xscale" ));
  }

} qfs_options_t, *qfs_options;

#ifdef  __cplusplus
extern "C" {
#endif

// Quasi-Scoring iteration
SEXP QSopt(SEXP R_start, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm, SEXP R_opt);

//SEXP invertMatrix(SEXP R_A, SEXP R_t);
//SEXP RSolve(SEXP R_A, SEXP R_B, SEXP R_type);

#ifdef  __cplusplus
}
#endif

#endif /* AUX_H_ */

