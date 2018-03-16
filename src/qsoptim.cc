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

#define STPMAX 10
#define TOLSTEP 1e-12

qfs_result qfscoring(double *x, int n, double &f, int &fntype,
			 qfs_options qfs, int &info);

void backtr(int n, double *xold, double &fold, double *d, double *g, double *x,
			double &f, int &check, int &fntype, double &stepmax, double &slope,
			 double &delta, double &rellen, qfs_options qfs, int &bounds, int pl, int &info);

// projection into box constraints
double med(double x, double y, double z, int &info);
void projmid(double *x, int nx, double *lb, double *ub, int &info);

// final status of QS iteration
SEXP getStatus(qfs_result status);

/** \brief
 *    Compute (projected) norm of quasi-score
 *    Components of QS are scaled and changed outside!
 *
 * @param x 	variable
 * @param data  data pointer
 * @param f 	pointer to function value
 * @param g 	pointer to gradient vector
 * @param info  integer, >0 if it reaches any bound constraints
 */
void fnQS(double *x, qfs_options qfs, double &f, int type, int &bounds, int &info) {
  ql_model qlm = qfs->qlm;
  int i=0, n=qlm->dx;

  projmid(x,n,qlm->lower,qlm->upper, bounds);
  if( (qlm->info = qlm->intern_qfScore(x)) != NO_ERROR) {		// compute score (always unscaled)
  	  WRR("Could not compute monitor function.")
  	 f=R_NaN;
  	 info=qlm->info;
  	 return;
  }
  if(type){
    f = qlm->intern_qfValue();								// no scaling for quasi-deviance
  } else {
	  double sum=0.,
			 *typf=qfs->typf,
			 *score=qlm->score;

	  for (i=0;i<n;++i)
	    sum += typf[i]*score[i]*typf[i]*score[i];    		// scaled norm^2 of quasi-score vector
	  if(!R_FINITE(sum)){
		  info=NaN_WARNING;
		  LOG_WARNING(info,"fnQS")
		  WRR("`NaN` detected in monitor function.")
	  }
	  f=0.5*sum;
  }
}

/*
 * Gradient of norm of quasi-score (fnQS)
 */
void fnGrad(qfs_options qfs, double *g, double *d, int type, int &info){
	ql_model qlm = qfs->qlm;
	int i=0, n=qlm->dx;
	double *typf=qfs->typf,					// scale quasi-score
		   *score=qlm->score;				// temporary

	if(type){
	  for(i=0;i<n;++i)
		d[i]=g[i]=score[i];						// simply use quasi-score as a gradient specification (no scaling)
	} else {
	  double *q=qlm->qlsolve.score;
	  for(i=0;i<n;++i){
	      d[i] = typf[i]*score[i];				// scaled d is used to compute Newton direction
	      q[i] = -typf[i]*d[i]; 				// gradient of norm^2 of quasi-score
	  }
	  matmult(qlm->qimat,n,n,q,n,1,g,info);
	  if(info > 0){
	    LOG_WARNING(info,"matmult")
 	    WRR("Could not compute gradient of monitor function.")
      }
	}
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
    	fntype=0, xdim = LENGTH(R_start);

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
    double fval=R_PosInf,
    	  qnorm=R_PosInf;

    qfs_result status = qfscoring(xsol,xdim,fval,fntype,&qfs,info);
    if(info)
     PRINT_MSG("Quasi-scoring iteration failed.")

    /* return objects */
    SEXP R_S, R_jac, R_I, R_Iobs=R_NilValue;
    PROTECT(R_S = allocVector(REALSXP,xdim));
    ++nProtected;
    PROTECT(R_I = allocMatrix(REALSXP,xdim,xdim));
    ++nProtected;
    PROTECT(R_jac = allocMatrix(REALSXP,xdim,qlm.nCov));
    ++nProtected;

    /* compute objectives */
    double *typf=qfs.typf,
    	   *score=qlm.score;

	if(fntype){
		double sum=0.0;
		for (int i=0;i<xdim;++i){
		  sum += SQR(typf[i]*score[i]);					/* scaling components of quasi-score vector */
		}
		if(!R_FINITE(sum)){
		  info=NaN_WARNING;
		  WRR("`NaN` detected in monitor function.")
		  LOG_WARNING(info,"Final computation failed of norm of quasi-score vector.")
		}
		qnorm = 0.5*sum;
	} else {
		qnorm = fval;									/* switch values: fval actually is norm^2 of quasi-score */
		fval = qlm.qfValue(qlm.score,qlm.qimat);
	}

    /* copy results  */
    MEMCPY(REAL(R_S),qlm.score,xdim);
    MEMCPY(REAL(R_I),qlm.qimat,xdim*xdim);
    MEMCPY(REAL(R_jac),qlm.jac,xdim*qlm.nCov);
    if(qfs.doIobs){
    	/* observed quasi-information matrix
    	 * as the Jacobian of the quasi-score vector (not scaled) */
      	PROTECT(R_Iobs = allocMatrix(REALSXP,xdim,xdim));
       	++nProtected;
       	info = qlm.intern_quasiObs(xsol, qlm.score, REAL(R_Iobs));
       	if( info != NO_ERROR )
       	  XWRR( info, "inter_quasiObs")
    }

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
    Rprintf("Qnorm: %f \n", qnorm);
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
     {"convergence", "message", "iter", "value", "par",
     "score", "sig2", "I", "Iobs", "varS", "start", "Qnorm",
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
    SET_VECTOR_ELT(R_ret, 8, R_Iobs);
    SET_VECTOR_ELT(R_ret, 9, R_varS);
    SET_VECTOR_ELT(R_ret, 10, R_start);
    SET_VECTOR_ELT(R_ret, 11, ScalarReal(qnorm));
    SET_VECTOR_ELT(R_ret, 12, mkString("qscoring"));
    SET_VECTOR_ELT(R_ret, 13, mkString("qle"));
    setVmatAttrib(&qlm, R_VmatNames, R_ret);
    SET_CLASS_NAME(R_ret,"QSResult")

    UNPROTECT(nProtected);
    return R_ret;
}


#define FREE_WORK { \
   FREE(g)     	    \
   FREE(d)			\
   FREE(xold)  	    \
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
					 double &f,  		 	/* objective value */
					 int &fntype,			/* =0 for norm^2 and =1 for quasi-deviance */
					 qfs_options qfs,    	/* options for scoring */
					 int &info)
{
   int i=0, niter=0, check=0, stopflag=0,
       bounds=0, Nmax=qfs->max_iter, pl=qfs->pl;

   // temp pointers
   ql_model qlm = qfs->qlm;
   double test=0.0, tmp=0.0, fold=0,
		  delta=1.0, slope=0.0, den=0.0, rellen=0.0,
		  *d=0, *xold=0, *g=0,
		  *typx=qfs->typx, *typf=qfs->typf,
		  *qimat=qlm->qimat, *score=qlm->score;

   CALLOCX(g,n,double);
   CALLOCX(d,n,double);
   CALLOCX(xold,n,double);

   /* first compute monitor */
   fnQS(x,qfs,f,fntype,bounds,info);
   if(info){
   	 FREE_WORK
   	 LOG_ERROR(NaN_ERROR,"`fnQS`: cannot evaluate monitor function.")
   	 return QFS_EVAL_ERROR;
   }
   /* and gradient (quasi-score vector if fntype > 0 */
   fnGrad(qfs,g,d,fntype,info);
   if(info){
	 FREE_WORK
	 LOG_ERROR(NaN_ERROR,"`fnGrad`: cannot evaluate gradient function.")
	 return QFS_EVAL_ERROR;
   }
   for (i=0;i<n; ++i) {
	 tmp=typf[i]*std::fabs(score[i]);
	 if (tmp > test) test=tmp;
   }
   if (test < 0.01*qfs->score_tol) {
	 FREE_WORK
     return QFS_SCORETOL_REACHED;
   }
   /* more strictive test first */
   if(f < 0.01*qfs->ftol_stop){
	FREE_WORK
    return QFS_STOPVAL_REACHED;
   }

   test=0.0;
   for (i=0;i<n;++i)
	 test += x[i]*typx[i]*x[i]*typx[i];
   double stepmax = STPMAX*MAX(std::sqrt(test),denorm(typx,n));
   if(!R_FINITE(stepmax)){
	   FREE_WORK
	   info=NaN_WARNING;
	   LOG_ERROR(info,"Cannot compute maximum step length")
	   return QFS_ERROR;
   }

   if(pl >= 10) {
	 Rprintf("stepmax.............%3.6f\n",stepmax);
   }

   fold = f;
   MEMCPY(xold,x,n);
   double *qlqimat = qlm->qlsolve.qimat;

   /*! optimization loop */
   for(niter=0; niter < Nmax; ++niter)
   {
	     /* solve for direction d = I^{-1} Q */
	     gsiSolve(qimat,n,d,1,qlqimat,info,Chol);
         if(info != 0){
        	 FREE_WORK
        	 PRINT_MSG("Cannot compute quasi-score correction vector.")
        	 LOG_ERROR(LAPACK_SOLVE_ERROR, "gsiSolve");
          	 qfs->num_iter=niter;
          	 return QFS_BAD_DIRECTION;
         }
         /* Line search:
          *  dynamically switch between both monitor functions */
         backtr(n,xold,fold,d,g,x,f,check,fntype,stepmax,slope,delta,rellen,qfs,bounds,pl,info);
         if(info > 0){
        	 FREE_WORK
			 qfs->num_iter=niter;
        	 LOG_ERROR(info,"backtr")
        	 return QFS_LINESEARCH_ERROR;
         }

         /*! display information */
         if(pl >= 10) {
           Rprintf("iteration......... %d \n", niter);
           Rprintf("objective......... %3.12f (fntype=%d) \n", f, fntype);
           Rprintf("at bounds......... %d \n", bounds);
           Rprintf("step max...........%3.4f \n", stepmax);
           Rprintf("step size..........%3.12f (check=%d) \n", delta, check);
           Rprintf("slope............ %3.12f \n\n", slope);
           printVector("par", x, &n);
           Rprintf("\n");
           printVector("score", score, &n);
           Rprintf("\n");
           printVector("direction", d, &n);
           Rprintf("length: %3.12f\n\n", denorm(d,n));
           printVector("gradient", g, &n);
           Rprintf("\n");
           Rprintf("--------------------------------------------------------------\n\n");
         }

         if(check > 0) {
        	 if(stopflag) {
        		FREE_WORK
				qfs->num_iter=niter;
        		if(check == 1)
        		 fnQS(x,qfs,f,fntype,bounds,info);      								/* re-compute objective  */
        		return (rellen < 1.5e-3 ? QFS_STEPMIN_REACHED : QFS_STEPTOL_REACHED);	/* rellen is length of scaled direction */
           	 } else {
        		fntype = (fntype > 0 ? 0 : 1);											/* type of monitor function */
        		fnQS(x,qfs,f,fntype,bounds,info);
        		if(info) break;
        		fnGrad(qfs,g,d,fntype,info);
        		if(info) break;
        		check=0;
        		stopflag=1;
        		continue;
        	 }
		 }
         /* update gradient at new x
          * and copy current score in d */
         fnGrad(qfs,g,d,fntype,info);
         if(info) break;

         /*! test for small enough objective */
         if(fntype) {
        	/*! test for score being zero */
			test=0.;
			for (i=0;i<n;++i) {
			 tmp = std::fabs(score[i]);
			 if(tmp > test) test=tmp;
			}
			if(test < qfs->score_tol) {
			  FREE_WORK
			  qfs->num_iter = niter;
			  return QFS_SCORETOL_REACHED;
			}
        	/* monitor is quasi-deviance */
            if(f < qfs->ftol_stop) {
              FREE_WORK
			  qfs->num_iter=niter;
              return QFS_STOPVAL_REACHED;
            }
            /* test for f relative change */
			if ( std::fabs(f-fold)/MAX(f,DBL_EPSILON) < qfs->ftol_rel) {
			  FREE_WORK
			  qfs->num_iter=niter;
			  return QFS_FTOLREL_REACHED;
			}

         } else {
        	/*! test for score being zero */
			test=0.;
			for (i=0;i<n;++i) {
			 tmp = typf[i] * std::fabs(score[i]);
			 if(tmp > test) test=tmp;
			}
			if(test < qfs->score_tol) {
			  FREE_WORK
			  qfs->num_iter = niter;
			  return QFS_SCORETOL_REACHED;
			}
        	/* monitor is norm^2 of quasi-score */
        	test=0.;
        	den=MAX(f,0.5*n);
        	/* relative gradient */
        	for (i=0;i<n;++i) {
				tmp=std::fabs(g[i])*MAX(std::fabs(x[i]),1./typx[i])/den;
				if(tmp > test) test=tmp;
			}
			/*! test for local min */
			if(test < qfs->grad_tol) {
				 FREE_WORK
				 qfs->num_iter=niter;
				 return (f < qfs->ftol_abs ? QFS_CONVERGENCE : QFS_LOCAL_CONVERGENCE);
			}
			/*! test for zero slope */
			if(std::fabs(slope) < qfs->slope_tol) {
				FREE_WORK
				qfs->num_iter=niter;
				return QFS_SLOPETOL_REACHED;
			}
         }
         /*! test for relative change in x */
		 test=0.;
		 for (i=0;i<n; ++i) {
		   tmp = (std::fabs(x[i]-xold[i])) / MAX(std::fabs(x[i]),1./typx[i]);
		   if(tmp > test) test = tmp;
		 }
		 if(test < qfs->xtol_rel) {
			FREE_WORK
			qfs->num_iter=niter;
			return QFS_XTOL_REACHED;
		 }
		 /* update */
		 fold=f;
		 MEMCPY(xold,x,n);
   } /*! end for */

   FREE_WORK
   qfs->num_iter=niter;
   return QFS_MAXITER_REACHED;
}

#undef FREE_WORK

void backtr(int n, double *xold, double &fold,  double *d, double *g, double *x,
		double &f, int &check,  int &fntype, double &stepmax, double &slope,
		 double &delta, double &rellen, qfs_options qfs, int &bounds, int pl, int &info) {

  int i=0;
  double s=1., tmp=0., dirlen=0.,
		 stepmin=0., alpha=1e-4, *typx=qfs->typx;

  info=0;
  check=0;
  for (i=0; i<n; ++i)
    tmp += typx[i] * d[i] * typx[i] * d[i];
  dirlen=std::sqrt(tmp);							/* newton length */
  if(!R_FINITE(dirlen)){
  	 LOG_ERROR(NaN_WARNING,"Cannot compute length of direction.")
  	 info=QFS_LINESEARCH_ERROR;
  	 return;
  }
  if (dirlen > stepmax){							/* not too big steps */
   for (i=0; i<n; ++i)
	d[i] *= (stepmax / dirlen);						/* scaled direction */
  }
  slope=(fntype ? -fold : innerProduct(g,d,n));     /* slope should be < 0 */
  if(!R_FINITE(slope) || slope >= 0.0){
	info=QFS_LINESEARCH_ZEROSLOPE;
	WRR("Roundoff error in line search.")
  	check=1;
  	return;
  }
  rellen=0.0;
  for (i=0; i<n; ++i) {
  	tmp=std::fabs(d[i])/MAX(std::fabs(xold[i]),1./typx[i]);
  	if (tmp > rellen) rellen = tmp;
  }
  if(rellen > 0.)
   stepmin=MAX(qfs->xtol_rel/rellen,TOLSTEP);
  else WRR("Relative step length should not be positive.")
  if(pl >= 100){
	  Rprintf("\n");
	  Rprintf("stepmin............%3.12f\n",stepmin);
	  Rprintf("rellength..........%3.12f\n\n",rellen);
	  printVector("scaled direction (d)",d,&n);
	  Rprintf("\n");
  }
  /* scaled direction length already too small */
  if(rellen < stepmin) {
	  check=2;
      return;
  }

  for(;;) {
	  delta=s;
	  for (i=0; i<n; ++i)
	  	x[i] = xold[i] + s * d[i];
	  fnQS(x,qfs,f,fntype,bounds,info);
	  if(info > 0) break;
	  if(f <= fold + alpha * s * slope) 			  /* set some required decrease */
	  	return;
	  if(s < stepmin){
		  f=fold;
		  for(i=0; i<n; ++i)
		    x[i] = xold[i];
		  check=1;
		  return;
	  } else s = .5*s;								/* decrease step */
  }
  /* error at this point */
  for (i=0;i<n;++i)
    x[i] = xold[i];
}

double med(double x, double y, double z, int &info) {
   if ( (x - y) * (z - x) >= 0 ) {
      if((x - y) * (z - x) == 0)
       info=1;
      return x;
   } else if ( (y - x) * (z - y) >= 0 ) {
       info=1;
       return y;
   } else {
       info=1;
       return z;
  }
}


void projmid(double *xproj, int nx, double *lb, double *ub, int &info) {
  int i=0;
  for(info=0;i<nx;++i)
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
       // (= +0)
       case QFS_CONVERGENCE:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_CONVERGENCE: Optimization stopped because of approximate root."));
           break;
       // (= +1)
       case QFS_SCORETOL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_SCORETOL_REACHED: Optimization stopped because score_tol was reached."));
           break;
       // (= +2)
       case QFS_FTOLREL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_FTOLREL_REACHED: Optimization stopped because ftol_rel was reached."));
           break;
       // (= +3)
       case QFS_STOPVAL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_STOPVAL_REACHED: Optimization stopped because ftol_stop  was reached."));
           break;
       // (= +4)
       case QFS_XTOL_REACHED:
           SET_STRING_ELT(R_message, 0, mkChar("QFS_XTOL_REACHED: Optimization stopped because xtol_rel was reached."));
           break;
   	   // (= +5)
	   case QFS_SLOPETOL_REACHED:
			SET_STRING_ELT(R_message, 0, mkChar("QFS_SLOPETOL_REACHED: Optimization stopped because slope_tol was reached."));
			break;
	   // (= +7)
	   case QFS_STEPMIN_REACHED:
	   	   	SET_STRING_ELT(R_message, 0, mkChar("QFS_STEPMIN_REACHED: Optimization stopped because minimum relative direction length was reached."));
	   	   	break;
	   // (= +8)
	   case QFS_STEPTOL_REACHED:
	        SET_STRING_ELT(R_message, 0, mkChar("QFS_STEPTOL_REACHED: Optimization stopped because minimum step length was reached."));
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
       case QFS_LINESEARCH_ERROR:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LINESEARCH_ERROR: Optimization stopped because sufficient decrease could not be found."));
            break;
       case QFS_LINESEARCH_ZEROSLOPE:
            SET_STRING_ELT(R_message, 0, mkChar("QFS_LINESEARCH_ZEROSLOPE: Optimization stopped because of nearly zero slope."));
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

/**
 *  @brief:
 *   testing internal function `invMatrix`,
 *   do not forget to reset dynamic symbols
 *
 *  */
//SEXP invertMatrix(SEXP R_A,SEXP R_type)
//{
//	int err = 0;
//	int n = GET_DIMS(R_A)[0];
//	inversion_type type = (inversion_type) asInteger(AS_INTEGER(R_type));
//
//	SEXP R_B = PROTECT(allocMatrix(REALSXP, n, n));
//	SEXP R_AP = PROTECT(allocMatrix(REALSXP, n, n));
//
//	MEMCPY(REAL(R_AP),REAL(R_A),SQR(n));
//	invMatrix(REAL(R_AP), n, REAL(R_B), &err, type);
//
//  UNPROTECT(2);
//	return R_B;
//}

/**
 * @brief
 * 	 testing internal function `gsiSolve`
 *   do not forget to reset dynamic symbols
 */
//SEXP RSolve(SEXP R_A, SEXP R_B, SEXP R_type){
//	int err = 0,
//	 n = GET_DIMS(R_A)[0],
//	 m = GET_DIMS(R_B)[1];
//
//	Rprintf("n: %d, m: %d \n",n,m);
//	Rprintf("A -> [%d,%d] \n",GET_DIMS(R_A)[0],GET_DIMS(R_A)[1]);
//	Rprintf("B -> [%d,%d] \n",GET_DIMS(R_B)[0],GET_DIMS(R_B)[1]);
//
//	inversion_type type = (inversion_type) asInteger(AS_INTEGER(R_type));
//
//	SEXP R_AP = PROTECT(allocMatrix(REALSXP, n, n));
//	MEMCPY(REAL(R_AP),REAL(R_A),n*n);
//
//	SEXP R_Z = PROTECT(allocMatrix(REALSXP, n, m));
//	MEMCPY(REAL(R_Z),REAL(R_B),n*m);
//
//	double *Awork = CALLOC(SQR(n),double);
//	gsiSolve(REAL(R_AP), n, REAL(R_Z), m, Awork, &err, type);
//	FREE(Awork);
//
//	UNPROTECT(2);
//	return R_Z;
//}

