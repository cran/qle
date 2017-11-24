/**
 * @file        optimize.cc
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions for optimization
 *
 * Explanation: Quasi-Fisher scoring iteration
 *
 */

#include "qsoptim.h"

#include <R_ext/Applic.h>
#include <R_ext/Constants.h>

#define TOLFAILS 200

static double ARG1, ARG2;

#define FMAX(a,b) (ARG1=(a),ARG2=(b),(ARG1) > (ARG2) ?\
        (ARG1) : (ARG2))

#define FMIN(a,b) (ARG1=(a),ARG2=(b),(ARG1) < (ARG2) ?\
        (ARG1) : (ARG2))

qfs_result qfscoring(double *x, int n, double *f,
			fnCompare fnMonitor, qfs_options cond, int *info);

void backtr(int n, double xold[], double fold, double p[], double x[],
			double *f, int *check, fnCompare monitor,
			double tau, double *delta, void *data);

double med(double x, double y, double z, int *info);

void projmid(double *x, int nx, double *lb, double *ub, int *info);

// projected objective function: score statistic
void fnLSProjMonitor(double *x, void *data, double *f);

SEXP getStatus( qfs_result status );

/** \brief Calculate Fisher score statistic
 *         of projected variable
 *
 * @param x variable
 * @param data data pointer
 * @param f pointer to function value
 */

void
fnLSProjMonitor(double *x, void *data, double *f) {
  qfs_options qfs = (qfs_options) data;
  projmid(x,qfs->qlm->glkm->dx,qfs->qlm->lower, qfs->qlm->upper, &(qfs->info));
  *f = qfs->qlm->intern_qfScoreStat(x);
}

/** \brief C Interface:
 *      Local root finding of score equations
 *
 * @param R_start start vector
 * @param R_args argument list for kriging
 * @param R_opt options for qfs iteration
 *
 * @return result object
 */

SEXP QSopt(SEXP R_start, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm, SEXP R_opt) {

    int nProtected = 0, info = 0,
		xdim = LENGTH(R_start);

    /*! set start vector */
    SEXP R_sol = R_NilValue;
    PROTECT(R_sol = allocVector(REALSXP, xdim));
    ++nProtected;

    double *xsol = REAL(R_sol);
    double *start = REAL(AS_NUMERIC(R_start));
    MEMCPY(xsol,start,xdim);

    ql_model_t qlm(R_qsd, R_qlopts, R_X, R_Vmat, R_cm, COPY_ZERO);

    /*! Set optimization options*/
    qfs_options_t qfs(&qlm,R_opt);

    /* scoring iteration */
    double fval = HUGE_VAL;
    fnCompare fnMonitor = &fnLSProjMonitor;
    qfs_result status = qfscoring(xsol,xdim,&fval,fnMonitor,&qfs,&info);

    /* return objects */
    SEXP R_S, R_jac, R_I;
    PROTECT(R_S = allocVector(REALSXP,xdim));
    ++nProtected;
    PROTECT(R_I = allocMatrix(REALSXP,xdim,xdim));
    ++nProtected;
    PROTECT(R_jac = allocMatrix(REALSXP,xdim,qlm.nCov));
    ++nProtected;

    /* copy results  */
    MEMCPY(REAL(R_S),qlm.score,xdim);
    MEMCPY(REAL(R_I),qlm.qimat,xdim*xdim);
    MEMCPY(REAL(R_jac),qlm.jac,xdim*qlm.nCov);

    /* add prediction variances to return list */
    SEXP R_sig2 = R_NilValue, R_varS = R_NilValue;

    // names variance matrix
    SEXP R_VmatNames;
    PROTECT(R_VmatNames = allocVector(VECSXP,2));
    ++nProtected;
    SEXP R_obs = getListElement( R_qsd, "obs" );
    SET_DIMNAMES_MATRIX(R_VmatNames,R_obs)

    if(qlm.glkm->krigType) {
      	PROTECT(R_sig2 = allocVector(REALSXP,qlm.nCov));
      	++nProtected;
      	MEMCPY(REAL(R_sig2),qlm.glkm->krigr[0]->sig2,qlm.nCov);
       	setAttrib(R_sig2,R_NamesSymbol,getAttrib(R_obs,R_NamesSymbol));
       	/* variance quasi-score */
       	PROTECT(R_varS = allocMatrix(REALSXP,qlm.dx,qlm.dx));
       	++nProtected;
       	qlm.intern_varScore(REAL(R_varS));
    }

#ifdef DEBUG
    Rprintf("value: %f \n", fval);
    printMatrix("vmat",qlm.qld->vmat, &xdim,&xdim);
    printMatrix("jac",REAL(R_jac), &xdim,&qlm.nCov);
    printMatrix("I",REAL(R_I), &xdim,&xdim);
    printVector("start:", xsol, &xdim);
    printVector("score:", REAL(R_S), &xdim);
#endif

    SEXP R_dimnames = R_NilValue;
    PROTECT(R_dimnames = allocVector(VECSXP,2));
    ++nProtected;
    SET_DIMNAMES_MATRIX(R_dimnames,R_start)
    setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
    setAttrib(R_sol, R_NamesSymbol, getAttrib(R_start,R_NamesSymbol));

    static const char *nms[] =
     {"status", "message", "iter", "value", "par",
     "score", "sig2", "I", "varS", "start", "convergence",
	 "method", "criterion", ""};

    SEXP R_ret = R_NilValue;
    PROTECT(R_ret = mkNamed(VECSXP, nms));
    ++nProtected;

    SET_VECTOR_ELT(R_ret, 0, ScalarInteger((int)status));
    SET_VECTOR_ELT(R_ret, 1, getStatus(status));
    SET_VECTOR_ELT(R_ret, 2, ScalarInteger(qfs.num_iter));
    SET_VECTOR_ELT(R_ret, 3, ScalarReal(fval));
    SET_VECTOR_ELT(R_ret, 4, R_sol);
    SET_VECTOR_ELT(R_ret, 5, R_S);
    SET_VECTOR_ELT(R_ret, 6, R_sig2);
    SET_VECTOR_ELT(R_ret, 7, R_I);
    SET_VECTOR_ELT(R_ret, 8, R_varS);
    SET_VECTOR_ELT(R_ret, 9, R_start);
    SET_VECTOR_ELT(R_ret, 10, ScalarInteger(info));
    SET_VECTOR_ELT(R_ret, 11, mkString("qscoring"));
    SET_VECTOR_ELT(R_ret, 12, mkString("qle"));
    setVmatAttrib(&qlm, R_VmatNames, R_ret);
    SET_CLASS_NAME(R_ret,"QSResult")

    UNPROTECT(nProtected);
    return R_ret;
}


#define FREE_WORK { \
   FREE(d)     	    \
   FREE(xold)  	    \
   FREE(gradf) 		\
}

/** \brief  Quasi-Fisher scoring iteration
 *        Comment: Either use step length equal to 1
 *                 or a line search based on the Fisher-Score statistic
 *
 * @param x start vector
 * @param n dimension of start vector
 * @param f monitor function value at solution
 * @param fnMonitor quasi-ddeviance function
 * @param qfs options
 * @param info information code
 *
 * @return result flag
 */

qfs_result qfscoring(double *x,			 	/* start */
					 int n,      		 	/* parameter length */
					 double *f,  		 	/* objective value */
					 fnCompare fnMonitor,	/* objective function */
					 qfs_options qfs,    	/* options for scoring */
					 int *info)
{

   /*! calculate value, score and qimat */
   fnMonitor(x, qfs, f);

   /*! test for f */
   if(*f < qfs->ftol_stop ){
     *info=QFS_CONVERGENCE;
     return QFS_STOPVAL_REACHED;
   }

   ql_model qlm = qfs->qlm;
   qfs_result status = QFS_NO_CONVERGENCE;

   int i=0, niter=0, check=0;
   int Nmax = qfs->max_iter, pl = qfs->pl;

   double fold=0 , test=0, tmp=0,
   	   	  slope2=0, tau=0.5, delta=1.0,
		  slope2_tol=SQR(qfs->slope_tol);

   double *d,*xold,*gradf;
   CALLOCX(d,n,double);
   CALLOCX(xold,n,double);
   CALLOCX(gradf,n,double);
   // score vector from `fnMonitor`
   double *score = qlm->score;

   test=0.0;
   for (i=0;i<n; ++i) {
      if (std::fabs(score[i]) > test)
        test=std::fabs(score[i]);
   }
   if (test < qfs->score_tol) {
      FREE_WORK
      *info=QFS_CONVERGENCE;
      return QFS_SCORETOL_REACHED;
    }

   /*! optimization loop */
   for(niter=0; niter < Nmax; ++niter)
   {
	     fold=*f;
	     MEMCPY(xold,x,n);
	     MEMCPY(gradf,qlm->score,n);
         solveCH(qlm->qimat,n,n,gradf,1,d,info);
         if(*info){
        	 WRR("Solving for the direction failed. Try something else.")
			 MEMCPY(d,gradf,n);
			 solveLU(qlm->qimat,n,d,n,info);
        	 XERR(*info,"Solving also failed by LU decomposition.");
         }

         backtr(n,xold,fold,d,x,f,&check,fnMonitor,tau,&delta,(void*) qfs);
         slope2 = innerProduct(d,d,n);

         /*! display information */
         if(pl >= 10) {
           Rprintf("--------------------------------------------------------------\n\n");
           Rprintf("iter:.............%d \n", niter);
           Rprintf("bounds............%d \n", qfs->info);
           Rprintf("value:............%3.10f \n", *f);
           Rprintf("LS step...........%3.10f (check=%d) \n\n", delta, check);
           printVector("par:", x, &n);
           Rprintf("\n");
           printVector("gradf:", gradf, &n);
           Rprintf("\n");
           printVector("score:", score, &n);
           Rprintf("\n");
           printVector("direction:", d, &n);
           Rprintf("\t norm_2: %3.10f\n\n", slope2);
         }

         /** test for zero slope of direction */
         if(slope2 < slope2_tol) {
             FREE_WORK
			 qfs->num_iter=niter;
             *info=(*f < qfs->ftol_abs ? QFS_CONVERGENCE : QFS_LOCAL_CONVERGENCE);
             status=QFS_SLOPETOL_REACHED;
             return status;
         }
         /*! test for score being zero */
         test=0.0;
         for (i=0;i<n;++i) {
        	tmp=std::fabs(score[i]);
            if(tmp > test)
              test=tmp;
         }
         if(test < qfs->score_tol) {
           FREE_WORK
		   qfs->num_iter = niter;
            *info=(*f < qfs->ftol_abs ? QFS_CONVERGENCE : QFS_LOCAL_CONVERGENCE);
           status=QFS_SCORETOL_REACHED;
           return status;
         }
         /*! test for f */
         if(*f < qfs->ftol_stop) {
             FREE_WORK
			 qfs->num_iter=niter;
             *info=QFS_CONVERGENCE;
             return QFS_STOPVAL_REACHED;
         }
         /*! test for local min */
         if(check>0) {
        	   test=0.0;
        	  *info=QFS_NO_CONVERGENCE;
			 status=QFS_LINESEARCH_FAILURE;

			 for (i=0;i<n;++i) {
				   tmp=std::fabs(gradf[i]);
				   if(tmp > test)
					  test=tmp;
			 }
			 // Rprintf("test gradient: %3.12f \n", test);
			 // Rprintf("\tnorm gradf: %3.12f\n", sqrt(innerProduct(gradf,gradf,n)));

			 if(test < qfs->grad_tol) {
				  *info=(*f < qfs->ftol_abs ? QFS_CONVERGENCE : QFS_LOCAL_CONVERGENCE);
				 status=QFS_GRADTOL_REACHED;
			 } else {
				MEMCPY(x,xold,n);
				fnMonitor(x, qfs, f);
			 }
			 FREE_WORK
			 qfs->num_iter=niter;
			 return status;
         } else {  /* line search success */
			 test=(std::fabs(*f-fold))/FMAX(std::fabs(*f),DBL_MIN);
			 if (test < qfs->ftol_rel) {
				 FREE_WORK
				 qfs->num_iter=niter;
				 *info=QFS_CONVERGENCE;
				 return QFS_FTOLREL_REACHED;
			 }
			 /*! test for relative change in x */
			 test=0.0;
			 for (i=0;i<n; ++i) {
				 tmp=(std::fabs(x[i]-xold[i]))/std::fabs(x[i]);
				 if(tmp > test)
				   test=tmp;
			 }
			 if(test < qfs->xtol_rel) {
				FREE_WORK
				qfs->num_iter=niter;
				*info = QFS_CONVERGENCE;
				return QFS_XTOL_REACHED;
			 }
         }
   } /*! end for */
   FREE_WORK
   *info=status;
   qfs->num_iter=niter;
   return QFS_MAXITER_REACHED;
}

#undef FREE_WORK

void backtr(int n, double *xold, double fold,  double *p, double *x, double *f, int *check,
             fnCompare monitor, double tau, double *delta, void *data) {
  int i=0, ntry=0;
  double s=1.0, stol=1e-10, EPS=1e-8;

  *check=0;
  for (i=0;i<n;++i)
    x[i]=xold[i]+p[i];
  monitor(x, data, f);

  while( (*f > fold - EPS) && s > stol && ntry < TOLFAILS ) {
       ++ntry;
       s=tau*s;
       for (i=0;i<n; ++i)
         x[i]=xold[i]+s*p[i];
       monitor(x, data, f);
  }
  if(*f > fold - EPS)
	*check=1;
  *delta=s;
}


double med(double x, double y, double z, int *info) {
   if ( (x - y) * (z - x) >= 0 ) {
      if((x - y) * (z - x) == 0)
       *info=1;
      return x;
   } else if ( (y - x) * (z - y) >= 0 ) {
       *info=1;
       return y;
   } else {
       *info=1;
       return z;
  }
}


void projmid(double *xproj, int nx, double *lb, double *ub, int *info) {
  *info=0;
  for(int i=0;i < nx; ++i)
    xproj[i] = med(xproj[i],lb[i],ub[i], info);
}

/** \brief Convert status message
 *
 * @param status
 * @return message
 */
SEXP getStatus( qfs_result status ) {
   // convert message to an R object
   SEXP R_message;
   PROTECT(R_message = allocVector(STRSXP, 1));

   switch ( status ) {
       // (= +1)
       case QFS_CONVERGENCE:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_CONVERGENCE: Optimization stopped because of approximate root."));
           break;
       // (= +2)
       case QFS_SCORETOL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_SCORETOL_REACHED: Optimization stopped because score_tol was reached."));
           break;
       // (= +3)
       case QFS_FTOLREL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_FTOLREL_REACHED: Optimization stopped because ftol_rel was reached."));
           break;
       // (= +4)
       case QFS_STOPVAL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_STOPVAL_REACHED: Optimization stopped because ftol_stop or ftol_abs was reached."));
           break;
       // (= +5)
       case QFS_XTOL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs was reached."));
           break;
       // (= +6)
       case QFS_GRADTOL_REACHED:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_GRAD_TOL_REACHED: Optimization stopped because grad_tol was reached."));
            break;
       // (= +7)
       case QFS_SLOPETOL_REACHED:
             SET_STRING_ELT(R_message, 0, mkChar("QFS_SLOPETOL_REACHED: Optimization stopped because slope_tol was reached."));
             break;
       // (= +10)
       case QFS_LOCAL_CONVERGENCE:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LOCAL_CONVERGENCE: Optimization stopped because of local convergence."));
            break;

       /* Error codes (negative return values): */

       // (= -1)
       case QFS_NO_CONVERGENCE:
		   SET_STRING_ELT(R_message, 0, mkChar("QFS_NO_CONVERGENCE: Optimization stopped because no convergence could be detected."));
		   break;
       // (= -2)
       case QFS_BAD_DIRECTION:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_FAILURE: Could not calculate search direction."));
           break;
       case QFS_LINESEARCH_FAILURE:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LINESEARCH_FAILURE: Optimization stopped because of line search failure."));
            break;
       case QFS_LINESEARCH_ZEROSLOPE:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LINESEARCH_ZEROSLOPE: Optimization stopped because of nearly zero slope during line search."));
            break;
       case QFS_ERROR:
		    SET_STRING_ELT(R_message, 0, mkChar("QFS_FAILURE: Generic failure code."));
		    break;
       case QFS_MAXITER_REACHED:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_MAXTIME_REACHED: Optimization stopped because maxiter (above) was reached."));
            break;
       default:
           SET_STRING_ELT(R_message, 0, mkChar("Unknown return status."));
           break;
       }

   UNPROTECT( 1 );
   return R_message;
}
