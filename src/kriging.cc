/**
 * @file        kriging.cc
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions for kriging estimation
 *
 * Explanation: Kriging statistics and derivatives
 *
 */

#ifndef KRIGING_CC_
#define KRIGING_CC_

#define FD_EPS 1e-4
#define KRIGE_TOLERANCE 1e-10

#include "kriging.h"

static int ONE_ELEMENT = 1;
static double ZERO_ELEMENT = 0;
static ql_model qlm_global = NULL ;

void intern_dualKmat(double *C, int Cdim, double *F, int Fdim, double *K);
//
void intern_dualKrigingWeights(double *Cmat, int nlx, double *Fmat, int fddim, double *data, double *w);
//
void solveKriging(double *Cinv,double *Fmat, double *X1mat, double *Qmat,
      int lx, int fddim, double *sigma2, double *w, double *sigma0, double *f0 );


/** \brief Internal function: Dual kriging matrix K,
 *              Chilès, J.-P. and Delfiner, P. (2008) Frontmatter, in Geostatistics:
 *              Modeling Spatial Uncertainty, John Wiley & Sons, Inc., Hoboken, NJ, USA.
 *              Chapter 3.4.8
 *
 * @param C covariance matrix
 * @param Cdim dimension of covariance matrix
 * @param F trend matrix
 * @param Fdim columns of F
 * @param K dual kriging matrix
 */
void
intern_dualKmat(double *C, int nc, double *F, int fddim, double *K) {
	int i=0, j=0, k=nc+fddim;

	for (j = 0; j < nc; j++)
	      for (i = 0; i < nc; i++)
			K[i+j*(nc+fddim)] = C[i+j*nc];

	for (j = 0; j < fddim; j++)
	      for (i = 0; i < nc; i++)
			K[nc*k+i+j*k] = F[i+j*nc];

	for (j = 0; j < fddim; j++)
	      for (i = 0; i < fddim; i++)
	                  K[nc*k+nc+i+j*k] = 0;

	for (j = 0; j < nc; j++)
		for (i = 0; i < fddim; i++)
        		K[(i+nc)+j*(nc+fddim)] = F[j+i*nc]; 	// F^t

}


/**\brief Interface to low level C functions
   *     - Get dual kriging weights for a single statistic
   *     - Comment: This function is independent of the
   *                underlying C data structure for the kriging models, which is not initialized
   *
   * @param R_Cmat Covariance matrix
   * @param R_Fmat Trend matrix
   * @param R_data Data frame of sampled statistic values at the design sites
   *
   * @return  R vector of kriging weights
   */

SEXP
getDualKrigingWeights(SEXP R_Cmat, SEXP R_Fmat, SEXP R_data) {
    int nx = GET_DIMS(R_Cmat)[0],
        fddim = GET_DIMS(R_Fmat)[1];
    int dK = nx+fddim;

    SEXP R_weights;
    PROTECT(R_weights = NEW_NUMERIC(dK));
    intern_dualKrigingWeights(REAL(R_Cmat), nx, REAL(R_Fmat), fddim, REAL(R_data), REAL(R_weights));

    UNPROTECT(1);
    return R_weights;

}


/**\brief Kriging prediction: calculate dual kriging weights
   *     - Kriging of general data
   *     - Comment: The kriging variance is not calculated
   *
   * @param Cmat Covariance matrix (with estimated parameters)
   * @param nlx Rows/Cols of Cmat
   * @param Fmat trend matrix
   * @param fddim Columns of trend matrix
   * @param data Data vector
   * @param w The weights
   *
   * return (void)
   */

void
intern_dualKrigingWeights(double *Cmat, int nlx, double *Fmat, int fddim, double *data, double *w) {
  int j=0, info=0, dK = nlx+fddim;

  double *K = 0;
  CALLOCX(K, dK*dK, double);
  intern_dualKmat(Cmat,nlx,Fmat,fddim,K);

  for (; j<dK; j++) w[j] = 0;
  MEMCPY(w,data,nlx);

  solveLU(K,dK,w,ONE_ELEMENT,&info);
  if(info != 0){
	  FREE(K);
	  XWRR(info,"LU decomposition/solving failed for dual kriging weights computation.")
  }
  FREE(K);
}

/**\brief Kriging prediction (dual)
 *
 *    Comment: Kriging variances not calculated
 *
 * @param s0           covariance vector between design sites and x0
 * @param f0           trend vector of x0
 * @param fddim        length of f0 (columns of trend matrix)
 * @param w  [IN]      kriging weights
 * @param zx [OUT]     predicted values: z*(x0)
 *
 * return (void)
 */

void
intern_dualKrigingPrediction(double *s0, double *f0, int fddim, double *w, int nz, double *mean) {
  int i = 0;
  double zx = 0;
  for(i = 0; i < nz; i++)
	zx += w[i]*s0[i];
  for(i = 0; i < fddim; i++)
	zx += w[i+nz]*f0[i];
  if(!R_finite(zx) || ISNA(zx) || ISNAN(zx))
	 WRR("`NaN` detected in computation of dual kriging prediction.")
  *mean = zx;
}

SEXP estimateJacobian(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_covT, SEXP R_krigtype) {
   if(!isMatrix(R_Xmat))
	 XERR(1,"Expected matrix of sample points.");

   int *dimX = GET_DIMS(R_Xmat),
		nCov = LENGTH(R_covT);

   int i = 0,
	   dx = dimX[1],
   	   npoints = LENGTH(R_points);

   glkrig_models glkm(R_covT,R_Xmat,R_data,R_krigtype,COPY_ZERO);

   SEXP R_names, R_jacobians, R_jac_element, R_jac;
   PROTECT(R_jacobians = allocVector(VECSXP,npoints));

   // names Jacobian
   SEXP R_dimT = R_NilValue;
   SEXP R_nI = VECTOR_ELT(R_points,0);
   PROTECT(R_dimT = allocVector(VECSXP,2));
   SET_DIMNAMES_MATRIX2(R_dimT,R_nI,R_data)

   double *point  = 0;
   //krig_result_s krigr(nCov,lx);
   for(; i < npoints; ++i) {
	  point = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
	  //printVector("point",point,&dx);
	  glkm.intern_kriging(point);

	  PROTECT(R_jac = allocMatrix(REALSXP, dx, nCov));
	  setAttrib(R_jac, R_DimNamesSymbol, R_dimT);

	  // get jacobian
	  glkm.intern_jacobian(point,REAL(R_jac));
	  PROTECT(R_jac_element = allocVector(VECSXP, 1));
	  PROTECT(R_names = allocVector(STRSXP, 1));
	  SET_STRING_ELT(R_names, 0, mkChar("mean"));
	  setAttrib(R_jac_element, R_NamesSymbol, R_names);

	  SET_VECTOR_ELT(R_jac_element,0,R_jac);
	  SET_VECTOR_ELT(R_jacobians,i,R_jac_element);
	  UNPROTECT(3);

   }
   UNPROTECT(2);
   return R_jacobians;
}


/**\brief C Interface: Kriging of multiple data
   *
   * @param R_Xmat      sample (design) matrix
   * @param R_data      list of observation vectors (as a data.frame)
   * @param R_points    matrix of unobserved points
   * @param R_covSList  list of covariance structures
   *
   * return List of points,
   *        each containes a list of kriging mean and variance
   *        and a list of kriging weights for each point
   */


SEXP kriging(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_covList, SEXP R_krigType)
{
     int i = 0, k = 0, npoints = 0,
         nCov = length(R_covList),
        *dimX = GET_DIMS(R_Xmat);

     int lx = dimX[0]; // rows, number of sample points

     /* is matrix or vector */
     if(isMatrix(R_points))
       npoints = GET_DIMS(R_points)[0];  // number of points to  predict gradient
     else
       ERR("expected matrix of sample points.")

     /* init all kriging models */
     glkrig_models glkm(R_covList,R_Xmat,R_data,R_krigType,COPY_ZERO);

     SEXP R_mean, R_sigma2, R_weights,
	 	  R_weights_tmp, R_tmp, R_retlist;

     PROTECT(R_retlist = allocVector(VECSXP, npoints));

     double *mean = 0, *sigma2 = 0,
    		*points = REAL(R_points);

     SEXP R_nT = getAttrib(R_data, R_NamesSymbol);

     if(glkm.krigType) {
    	 const char *nms[] = {"mean", "sigma2", "weights", ""};
    	 for(; i<npoints; i++, points++) {
             PROTECT(R_mean = allocVector(REALSXP, nCov));
        	 PROTECT(R_sigma2 = allocVector(REALSXP, nCov));
             PROTECT(R_weights_tmp = allocVector(VECSXP, nCov));

             mean = REAL(R_mean);
             sigma2 = REAL(R_sigma2);
             for(k=0; k < nCov; ++k)   {
                 PROTECT(R_weights = allocVector(REALSXP, lx));
                 glkm.km[k]->univarKriging(points,npoints,mean+k,sigma2+k,REAL(R_weights));
                 SET_VECTOR_ELT(R_weights_tmp, k, R_weights);
                 UNPROTECT(1);
             }
             setAttrib(R_mean,R_NamesSymbol,R_nT);
             setAttrib(R_sigma2,R_NamesSymbol,R_nT);

             PROTECT(R_tmp = mkNamed(VECSXP, nms));
             SET_VECTOR_ELT(R_tmp, 0, R_mean);
             SET_VECTOR_ELT(R_tmp, 1, R_sigma2);
             SET_VECTOR_ELT(R_tmp, 2, R_weights_tmp);

             SET_VECTOR_ELT(R_retlist, i, R_tmp);
             UNPROTECT(4);
    	 }
      } else {
    	  const char *nms[] = {"mean", ""};
    	  for(; i<npoints; i++, points++) {
    		 PROTECT(R_mean = allocVector(REALSXP, nCov));

    		 mean = REAL(R_mean);
             for(k=0; k < nCov; k++)
               glkm.km[k]->dualKriging(points,npoints,mean+k);

             setAttrib(R_mean,R_NamesSymbol,R_nT);
             PROTECT(R_tmp = mkNamed(VECSXP, nms));
             SET_VECTOR_ELT(R_tmp, 0, R_mean);
             SET_VECTOR_ELT(R_retlist, i, R_tmp);
             UNPROTECT(2);
    	  }
     }

     SET_CLASS_NAME(R_retlist,"krigResult");

     UNPROTECT(1);
     return R_retlist;
}


/** \brief Solve the kriging equations, see
 *     Chilès, J.-P. and Delfiner, P. (2008) Frontmatter, in Geostatistics:
 *              Modeling Spatial Uncertainty, John Wiley & Sons, Inc., Hoboken, NJ, USA.
 *              page: 169
 *
 * @param Cinv  inverse covariance matrix (stored)
 * @param Fmat  trend matrix
 * @param X1mat X1 matrix, see function univarKriging
 * @param Qmat  Q matrix,  see function univarKriging
 * @param Nx rows of sample matrix
 * @param fddim columns of trend matrix
 * @param sigma2 prediction (kriging) variance
 * @param w kriging weights
 * @param sigma0 covariance vector between x0 (prediction point) and sample locations
 * @param f0 trend vector at x0 prediction point
 */

void solveKriging(double *Cinv,double *Fmat, double *X1mat, double *Qmat,
		 int lx, int fddim, double *sigma2, double *lambda,
		 double *sigma0, double *f0 ) {

	/* solving Kriging equations
	 *    lambdak = C^-1 sigma0
	 *    R = F^t lambdak -f0
	 *    mu = Q^-1 R
	 *    lambda = lambdak - X1 mu
	 */

    int info = 0, j = 0, k = 0;

	double *lambdak = 0,
			*Rvec = 0,
			*mu = 0;

	CALLOCX(lambdak,lx, double);
    CALLOCX(Rvec,fddim,double);
	CALLOCX(mu,fddim, double);

	matmult(Cinv,lx,lx,sigma0,lx,ONE_ELEMENT,lambdak,&info);
	if(info)
	  WRR("Na detected in matrix multiplication")
	matmult_trans(Fmat,lx,fddim,lambdak,lx,ONE_ELEMENT,Rvec,&info);
	if(info)
	  WRR("Na detected in matrix multiplication")

	for(k=0; k<fddim; k++)
	  mu[k] = Rvec[k] - f0[k];

	/*!  solve Q mu = R */
	solve_DSPTRS(Qmat,fddim,mu,ONE_ELEMENT,&info);
	if(info != 0){
	  FREE(lambdak);
	  FREE(Rvec);
	  FREE(mu);
	  *sigma2 = R_NaN;
	  *lambda = R_NaN;
	  XWRR(info,"Solving in kriging procedure failed.");
	  return;
	}
	matmult(X1mat, lx, fddim, mu, fddim, ONE_ELEMENT, lambda,&info);
	//Rprintf("%s:%u: %f", __FILE__, __LINE__, *sigma2);

	for(j=0; j < lx; j++)
	  lambda[j] = lambdak[j] - lambda[j];

	for (j=0; j < lx; j++)
	  *sigma2 -= lambda[j] * sigma0[j];

	for (j=0; j < fddim; j++)
	  *sigma2 -= mu[j] * f0[j];

	if (*sigma2<0) {
	    /* for numerical stability */
	    *sigma2 = KRIGE_TOLERANCE;
	} else if(*sigma2<KRIGE_TOLERANCE) {
	    *sigma2 = 0.0;
	}
	//Rprintf("%s:%u: %f", __FILE__, __LINE__, *sigma2);
	FREE(lambdak);
	FREE(Rvec);
	FREE(mu);
}

SEXP initQL(SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm)
{
  if(qlm_global) {
	 WRR("qldata is already initialized -> finalizing.");
	 if(!finalizeQL())
	   ERR("Could not free memory of `qlm_model`.")
  }
  if( (qlm_global = new(std::nothrow) ql_model_s(R_qsd, R_qlopts, R_X, R_Vmat, R_cm, COPY_ONE)) == NULL)
  	MEM_ERR(1,ql_model_s);

  return ScalarLogical(TRUE);
}



/**
 * @brief  FREE 'qldata' storage from R level
 * @return NULL
 */
// TODO:
SEXP finalizeQL() {
  if(qlm_global)
	DELETE(qlm_global)
  return ScalarLogical(TRUE);
}


SEXP qDValue(SEXP R_point) {
  if(!qlm_global) {
     ERR("Pointer to `qldata` object not set (NULL).");
     return R_NilValue;
  }
  return ScalarReal(qlm_global->intern_qfScoreStat(REAL(AS_NUMERIC(R_point))));
}

/**
 * \brief
 * 		Calculate Mahalanobis distance of statistics.
 *
 * 		If R_Sigma is used as inverse covariance matrix of statistics
 * 		(simulated), then the score vector as the gradient is returned,
 * 		otherwise the covariance matrix is continuously computed at each point acc.
 * 		to the settings of `ql`.
 */
SEXP mahalValue(SEXP R_point) {
  if(!qlm_global) {
     ERR("Pointer to `qldata` object not set (NULL).");
     return R_NilValue;
  }
  int info = 0,
	  dx = qlm_global->glkm->dx,
	  nCov = qlm_global->glkm->nCov;

  double f = 0;
  double *x = REAL(AS_NUMERIC(R_point));

  if(qlm_global->qld->qlopts.useSigma) {
  	 f = qlm_global->intern_mahalValue(x);
  } else {
     f = qlm_global->intern_mahalValue_theta(x);
  }

  SEXP R_ans, R_score;
  PROTECT(R_ans = ScalarReal(f));
  PROTECT(R_score = allocVector(REALSXP,dx));

  /* Jacobian */
  qlm_global->glkm->intern_jacobian(x,qlm_global->jac);

  /* score vector */
  matmult(qlm_global->jac,dx,nCov,qlm_global->qld->tmp,nCov,ONE_ELEMENT,REAL(R_score),&info);
  setAttrib(R_ans,install("score"),R_score);

  UNPROTECT(2);
  return R_ans;
}

/* add variance matrix to result list as
 * an attribute if kriging is used for
 * variance matrix approximation
*/
void setVmatAttrib(ql_model qlm, SEXP R_VmatNames, SEXP R_ans) {
    SEXP R_Vmat = R_NilValue;
    PROTECT(R_Vmat = allocMatrix(REALSXP,qlm->nCov,qlm->nCov));
    MEMCPY(REAL(R_Vmat),qlm->qld->vmat,qlm->nCov2);

    setAttrib(R_Vmat, R_DimNamesSymbol, R_VmatNames);
    setAttrib(R_ans, install("Sigma"), R_Vmat);
    UNPROTECT(1);
}

/**
 * \brief
 *
 * 	  Mahalanobis distance
 *
 *    1. Using a constant (already inverted) Sigma
 *    2. Average approximations: with and without prediction variances
 *    3. Kriging approximation of variance matrix
 *
 *    either returns only the value or with computed components
 */

SEXP mahalanobis(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm, SEXP R_qdValue)
{
	  int i = 0, info = 0,
		 np = LENGTH(R_points);
	  value_type type = (value_type) asInteger(AS_INTEGER(R_qdValue));

	  ql_model_t qlm(R_qsd, R_qlopts, R_X, R_Vmat, R_cm, type);

	  GLKM glkm = qlm.glkm;
	  ql_data qld = qlm.qld;
	  ql_options opts = &qld->qlopts;

	  int dx = glkm->dx,
		  nCov = glkm->nCov;

      SEXP R_ret, R_ans,
	  	   R_S, R_jac, R_I;
	  SEXP R_varS = R_NilValue,
		   R_sig2 = R_NilValue;

      /* constant Sigma
       *
       * Comment:
       *
       * Constant variance matrix is stored
       * as an attribute at the corresponding R function
       * */
	  if(opts->useSigma)
	  {
		  /* only values */
		  if(type) {
			  PROTECT(R_ret = allocVector(REALSXP,np));
			  double *fx = REAL(R_ret);

			  if(type == COPY_TRACE) {
				  for(i=0; i < np; ++i)
					fx[i] = qlm.intern_mahalVarTrace(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
			  } else {
 				  for(i=0; i < np; ++i)
			        fx[i] = qlm.intern_mahalValue(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
			  }

		  } else {

			  PROTECT(R_ret = allocVector(VECSXP,np));
			  // names QI
			  SEXP R_nI = VECTOR_ELT(R_points,0);
			  SEXP R_dimnames = R_NilValue;
			  PROTECT(R_dimnames = allocVector(VECSXP,2));
			  SET_DIMNAMES_MATRIX(R_dimnames,R_nI)
			  // names Jacobian
			  SEXP R_dimT = R_NilValue;
			  SEXP R_nT = getListElement( R_qsd, "obs" );
			  PROTECT(R_dimT = allocVector(VECSXP,2));
			  SET_DIMNAMES_MATRIX2(R_dimT,R_nI,R_nT)

			  /* using prediction variances */
			  if(glkm->krigType)
			  {
            	  double fval = 0;
            	  const char *nms[] = {"value", "par", "I", "score", "sig2", "jac","varS",""};

				  for(i=0;i<np;i++)
				  {
					 PROTECT(R_S = allocVector(REALSXP,dx));
					 PROTECT(R_sig2 = allocVector(REALSXP,nCov));
					 PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
					 PROTECT(R_varS = allocMatrix(REALSXP,dx,dx));
					 PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));

					 double *x = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));

					 /* mahalanobis distance */
					 fval = qlm.intern_mahalValue(x);
					 glkm->intern_jacobian(x,REAL(R_jac));

					 /* score vector */
					 matmult(REAL(R_jac),dx,nCov,qlm.qld->tmp,nCov,ONE_ELEMENT,REAL(R_S),&info);

					 /* quasi-info */
					 mat_trans(qld->jactmp,nCov,REAL(R_jac),dx,dx,nCov);
					 matmult(qld->vmat,nCov,nCov,qld->jactmp,nCov,dx,qld->Atmp,&info);
					 matmult(REAL(R_jac),dx,nCov,qld->Atmp,nCov,dx,REAL(R_I),&info);

					 /* We need prediction variances
					  * for the variance of quasi-score. The variance
					  * matrix is fixed and thus has no additional
					  * diagonal terms, i.e. kriging/CV variances)
					  */
					 qlm.intern_cvError(x);
					 qlm.intern_varScore(REAL(R_varS));
					 /* copy prediction variance */
					 MEMCPY(REAL(R_sig2),qlm.glkm->krigr[0]->sig2,nCov);

					 /*  set names but not for varS! */
					 setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
					 setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
					 setAttrib(R_sig2,R_NamesSymbol,getAttrib(R_nT,R_NamesSymbol));

					 PROTECT(R_ans = mkNamed(VECSXP, nms));
					 SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
					 SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
					 SET_VECTOR_ELT(R_ans, 2, R_I);
					 SET_VECTOR_ELT(R_ans, 3, R_S);
					 SET_VECTOR_ELT(R_ans, 4, R_sig2);
					 SET_VECTOR_ELT(R_ans, 5, R_jac);
					 SET_VECTOR_ELT(R_ans, 6, R_varS);

					 SET_VECTOR_ELT(R_ret, i, R_ans);
					 UNPROTECT(6);
				  }

			  } else {

				  /* no prediction variances */
				  double fval=0;
				  const char *nms[] = {"value", "par", "I", "score", "jac",""};

				  for(i=0;i<np;i++)
				  {
					 PROTECT(R_S = allocVector(REALSXP,dx));
					 PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
					 PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));

					 double *x = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
					 /* mahalanobis distance */
					 fval = qlm.intern_mahalValue(x);
					 glkm->intern_jacobian(x,REAL(R_jac));

					 /* score vector */
					 matmult(REAL(R_jac),dx,nCov,qld->tmp,nCov,ONE_ELEMENT,REAL(R_S),&info);
					 /* quasi-info */
					 mat_trans(qld->jactmp,nCov,REAL(R_jac),dx,dx,nCov);
					 matmult(qld->vmat,nCov,nCov,qld->jactmp,nCov,dx,qld->Atmp,&info);
					 matmult(REAL(R_jac),dx,nCov,qld->Atmp,nCov,dx,REAL(R_I),&info);

					 setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
					 setAttrib(R_I, R_DimNamesSymbol, R_dimnames);

					 PROTECT(R_ans = mkNamed(VECSXP, nms));
					 SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
					 SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
					 SET_VECTOR_ELT(R_ans, 2, R_I);
					 SET_VECTOR_ELT(R_ans, 3, R_S);
					 SET_VECTOR_ELT(R_ans, 4, R_jac);

					 SET_VECTOR_ELT(R_ret, i, R_ans);
					 UNPROTECT(4);
				  }
			  }
			  UNPROTECT(2); // matrix dimnames
		  }

	  } else {

		  if(type) {
			  PROTECT(R_ret = allocVector(REALSXP,np));
  			  double *f = REAL(R_ret);

  			  if(type == COPY_TRACE) {
			    for(i=0; i < np; ++i)
				  f[i] = qlm.intern_qfTrace(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
  			  } else {
  				for(i=0;i<np;i++)
  				  f[i] = qlm.intern_mahalValue_theta(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
  			  }

		  } else {

			  PROTECT(R_ret = allocVector(VECSXP,np));
			  // names QI
			  SEXP R_nI = VECTOR_ELT(R_points,0);
			  SEXP R_dimnames = R_NilValue;
			  PROTECT(R_dimnames = allocVector(VECSXP,2));
			  SET_DIMNAMES_MATRIX(R_dimnames,R_nI)
			  // names Jacobian
			  SEXP R_dimT = R_NilValue;
			  SEXP R_nT = getListElement( R_qsd, "obs" );
			  PROTECT(R_dimT = allocVector(VECSXP,2));
			  SET_DIMNAMES_MATRIX2(R_dimT,R_nI,R_nT)
			  // names variance matrix
			  SEXP R_VmatNames = R_NilValue;
			  PROTECT(R_VmatNames = allocVector(VECSXP,2));
			  SET_DIMNAMES_MATRIX(R_VmatNames,R_nT)

			  if(glkm->krigType)
			  {
					  double fval=0;
					  const char *nms[] = {"value", "par", "I", "score", "sig2", "jac","varS",""};

					  for(i=0;i<np;i++)
					  {
							R_CheckUserInterrupt();
							PROTECT(R_S = allocVector(REALSXP,dx));
							PROTECT(R_sig2 = allocVector(REALSXP,nCov));
							PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
							PROTECT(R_varS = allocMatrix(REALSXP,dx,dx));
							PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));

							double *x = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
							fval = qlm.intern_mahalValue_theta(x);

							// Jacobian
							glkm->intern_jacobian(x,REAL(R_jac));
							// score vector
							matmult(REAL(R_jac),dx,nCov,qld->tmp,nCov,ONE_ELEMENT,REAL(R_S),&info);

							// quasi-Info
							qlm.intern_quasiInfo(REAL(R_jac),REAL(R_I));

							//* variance score vector */
							qlm.intern_varScore(REAL(R_varS));
							/* copy prediction variance */
							MEMCPY(REAL(R_sig2),qlm.glkm->krigr[0]->sig2,nCov);
							/* set names */
							setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
							setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
							setAttrib(R_sig2,R_NamesSymbol,getAttrib(R_nT,R_NamesSymbol));

							PROTECT(R_ans = mkNamed(VECSXP, nms));
							SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
							SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
							SET_VECTOR_ELT(R_ans, 2, R_I);
							SET_VECTOR_ELT(R_ans, 3, R_S);
							SET_VECTOR_ELT(R_ans, 4, R_sig2);
							SET_VECTOR_ELT(R_ans, 5, R_jac);
							SET_VECTOR_ELT(R_ans, 6, R_varS);
							setVmatAttrib(&qlm, R_VmatNames, R_ans);

							SET_VECTOR_ELT(R_ret, i, R_ans);
							UNPROTECT(6);
					  }

			  } else {

					  double fval=0;
					  const char *nms[] = {"value", "par", "I", "score", "jac",""};

					  for(i=0;i<np;i++)
					  {
							R_CheckUserInterrupt();
							PROTECT(R_S = allocVector(REALSXP,dx));
							PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
							PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));

							double *x = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
							fval = qlm.intern_mahalValue_theta(x);

							// Jacobian
							glkm->intern_jacobian(x,REAL(R_jac));
							// score vector
							matmult(REAL(R_jac),dx,nCov,qld->tmp,nCov,ONE_ELEMENT,REAL(R_S),&info);

							// quasi-Info
							qlm.intern_quasiInfo(REAL(R_jac),REAL(R_I));

							setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
							setAttrib(R_I, R_DimNamesSymbol, R_dimnames);

							PROTECT(R_ans = mkNamed(VECSXP, nms));
							SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
							SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
							SET_VECTOR_ELT(R_ans, 2, R_I);
							SET_VECTOR_ELT(R_ans, 3, R_S);
							SET_VECTOR_ELT(R_ans, 4, R_jac);
							setVmatAttrib(&qlm, R_VmatNames, R_ans);

							SET_VECTOR_ELT(R_ret, i, R_ans);
							UNPROTECT(4);
					  }
			  }
			  UNPROTECT(3); // matrix dimnames
		  }
	}

	UNPROTECT(1);
	return R_ret;

}

/* inverse variance matrix has
 *  to be computed for each theta!
 *
 *  Only call if not using Sigma!
 *
 */
double
ql_model_s::intern_mahalValue_theta(double *x) {
  int k = 0;
  /* kriging */
  glkm->intern_kriging(x);
  krig_result krig = glkm->krigr[0];

  /* continuously update variance matrix */
  /* without calculation of score vector as the gradient */
  for(k = 0; k < nCov; ++k)
	qld->tmp[k] = qld->qtheta[k] = qld->obs[k] - krig->mean[k];

  intern_cvError(x);
  varMatrix(x,krig->sig2,qld->vmat);

  solve_DSPTRS(qld->vmat,nCov,qld->tmp,ONE_ELEMENT, &info);
  if(info != 0){
	XWRR(info,"Internal computation of Mahalanobis distance failed.")
	return R_NaN;
  }
  double sum = 0;
  for(k = 0; k < nCov; ++k)
	 sum += qld->tmp[k] * qld->qtheta[k];
  return sum;
}

/* vmat is inverse variance matrix */
double
ql_model_s::intern_mahalValue(double *x) {
	int k = 0, info = 0;
	glkm->intern_kriging(x);
	krig_result krig = glkm->krigr[0];
 	for(k = 0; k < nCov; ++k)
	   qld->qtheta[k] = qld->obs[k] - krig->mean[k];

	matmult(qld->vmat,nCov,nCov,qld->qtheta,nCov,ONE_ELEMENT,qld->tmp,&info);
	if(info != 0){
	  XWRR(info,"Internal computation of Mahalanobis distance failed.")
	  return R_NaN;
	}
	double sum=0;
	for(k = 0; k < nCov; ++k)
	  sum += qld->tmp[k] * qld->qtheta[k];
	return sum;
}

double
ql_model_s::intern_mahalVarTrace(double *x) {
	 int k = 0, info = 0;
	 double sum = 0;

	 glkm->intern_kriging(x);
	 krig_result krig = glkm->krigr[0];
	 for(k = 0; k < nCov; ++k)
	   qld->qtheta[k] = qld->obs[k] - krig->mean[k];

	 matmult(qld->vmat,nCov,nCov,qld->qtheta,nCov,ONE_ELEMENT,qld->tmp,&info);

	 /* score vector */
	 glkm->intern_jacobian(x,jac);
	 matmult(jac,dx,nCov,qld->tmp,nCov,ONE_ELEMENT,score,&info);
	 if(info != 0){
	 	   XWRR(info,"NaN detected in matrix multiplication.")
	 	   return R_NaN;
	 }
	 /* variance of score vector stored in qimat */
	 mat_trans(qld->jactmp,nCov,jac,dx,dx,nCov);
	 matmult(qld->vmat,nCov,nCov,qld->jactmp,nCov,dx,qld->Atmp,&info);
	 if(info != 0){
	   XWRR(info,"NaN detected in matrix multiplication.")
	   return R_NaN;
	 }
	 /* We need prediction variances
	  * for the variance of quasi-score. The variance
	  * matrix is fixed and thus has no additional
	  * diagonal terms, i.e. kriging/CV variances)
	  */
	 intern_cvError(x);
	 intern_varScore(qimat);

	 for(k = 0; k < dx; ++k)
	     sum += qimat[k*dx+k];
	 return sum/(double)dx;
}

/**
 *  \brief A test implementation for
 *  	   Cross-validation based prediction errors
 *
 */
//SEXP cvError(SEXP R_point, SEXP R_Xmat, SEXP R_data, SEXP R_covT, SEXP R_krigType, SEXP R_cm)
//{
//	/* CV models (might have different number of left out points! */
//	if(!isMatrix(R_Xmat))
//	 ERR("Expected sample matrix.");
//
//	cv_model_s cvmod(R_Xmat,R_data,R_cm);
//	/* full kriging models */
//	glkrig_models glkm(R_covT,R_Xmat,R_data,R_krigType,FALSE);
//
//	SEXP R_ret = PROTECT(allocVector(REALSXP,glkm.nCov));
//	cvmod.cvError(REAL(AS_NUMERIC(R_point)),glkm.km,REAL(R_ret));
//
//	UNPROTECT(1);
//	return R_ret;
//}


SEXP quasiDeviance(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm, SEXP R_qdValue)
{
      int i = 0, np = LENGTH(R_points);
      value_type type = (value_type)asInteger(AS_INTEGER(R_qdValue));

      ql_model_t qlm(R_qsd, R_qlopts, R_X, R_Vmat, R_cm, type);

      if(type) {
    	  SEXP Rval;
    	  PROTECT(Rval = allocVector(REALSXP,np));
          double *fx = REAL(Rval);

          if(type > COPY_ONE && qlm.glkm->krigType) {
        	 if(type == COPY_MOD){
				  for(;i < np; ++i)
				   fx[i] = qlm.intern_qfVarStat(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
        	 } else {
			   for(;i < np; ++i)
				 fx[i] = qlm.intern_qfTrace(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
			}
  		  } else {
  			for(; i < np; ++i)
  			  fx[i] = qlm.intern_qfScoreStat(REAL(AS_NUMERIC(VECTOR_ELT(R_points,i))));
		  }
          UNPROTECT(1);
          return Rval;

      } else {

    	  int dx = qlm.dx,
    		  nCov = qlm.nCov;

    	  double *point = NULL,
    			  fval = 0, qval = 0;

          SEXP R_ret, R_ans, R_sig2 = R_NilValue;
          SEXP R_Iobs, R_S, R_jac, R_I, R_varS;
          PROTECT(R_ret = allocVector(VECSXP,np));

          // names QI
          SEXP R_nI = VECTOR_ELT(R_points,0);
          SEXP R_dimnames = R_NilValue;
          PROTECT(R_dimnames = allocVector(VECSXP,2));
          SET_DIMNAMES_MATRIX(R_dimnames,R_nI)

          // names Jacobian
          SEXP R_dimT = R_NilValue;
          SEXP R_nT = getListElement( R_qsd, "obs" );
          PROTECT(R_dimT = allocVector(VECSXP,2));
          SET_DIMNAMES_MATRIX2(R_dimT,R_nI,R_nT)

          // names variance matrix
          SEXP R_VmatNames = R_NilValue;
          PROTECT(R_VmatNames = allocVector(VECSXP,2));
          SET_DIMNAMES_MATRIX(R_VmatNames,R_nT)

          if(qlm.glkm->krigType) {
        	  const char *nms[] =
        	  	  {"value", "par", "I", "score", "sig2",
        	  	   "jac","varS", "Iobs", "qval", ""};

        	  for(; i < np; ++i) {
					  PROTECT(R_S = allocVector(REALSXP,dx));
					  PROTECT(R_sig2 = allocVector(REALSXP,nCov));
					  PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
					  PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));
					  PROTECT(R_varS = allocMatrix(REALSXP,dx,dx));
					  PROTECT(R_Iobs = allocMatrix(REALSXP,dx,dx));

					  point = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
					  fval = qlm.qfScoreStat(point,REAL(R_jac),REAL(R_S),REAL(R_I));
					  qlm.intern_varScore(REAL(R_varS));
					  qlm.intern_quasiObs(point,REAL(R_S),REAL(R_Iobs));
					  qval = qlm.qfValue(REAL(R_S),REAL(R_varS));
					  MEMCPY(REAL(R_sig2),qlm.glkm->krigr[0]->sig2,nCov);

					  setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
					  setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
					  setAttrib(R_Iobs, R_DimNamesSymbol, R_dimnames);
					  setAttrib(R_sig2,R_NamesSymbol,getAttrib(R_nT,R_NamesSymbol));

					  PROTECT(R_ans = mkNamed(VECSXP, nms));
					  SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
					  SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
					  SET_VECTOR_ELT(R_ans, 2, R_I);
					  SET_VECTOR_ELT(R_ans, 3, R_S);
					  SET_VECTOR_ELT(R_ans, 4, R_sig2);
					  SET_VECTOR_ELT(R_ans, 5, R_jac);
					  SET_VECTOR_ELT(R_ans, 6, R_varS);
					  SET_VECTOR_ELT(R_ans, 7, R_Iobs);
					  SET_VECTOR_ELT(R_ans, 8, ScalarReal(qval));
					  setVmatAttrib(&qlm, R_VmatNames, R_ans);

					  SET_VECTOR_ELT(R_ret, i, R_ans);
					  UNPROTECT(7);
        	  }
          } else {
        	  const char *nms[] =
				  {"value", "par", "I", "score", "jac", "Iobs" ,""};

			  for(; i < np;  ++i) {
					  PROTECT(R_S = allocVector(REALSXP,dx));
					  PROTECT(R_I = allocMatrix(REALSXP,dx,dx));
					  PROTECT(R_jac = allocMatrix(REALSXP,dx,nCov));
					  PROTECT(R_Iobs = allocMatrix(REALSXP,dx,dx));

					  point = REAL(AS_NUMERIC(VECTOR_ELT(R_points,i)));
					  fval = qlm.qfScoreStat(point,REAL(R_jac),REAL(R_S),REAL(R_I));
					  qlm.intern_quasiObs(point,REAL(R_S),REAL(R_Iobs));

					  setAttrib(R_jac, R_DimNamesSymbol, R_dimT);
					  setAttrib(R_I, R_DimNamesSymbol, R_dimnames);
					  setAttrib(R_Iobs, R_DimNamesSymbol, R_dimnames);

					  PROTECT(R_ans = mkNamed(VECSXP, nms));
					  SET_VECTOR_ELT(R_ans, 0, ScalarReal(fval));
					  SET_VECTOR_ELT(R_ans, 1, VECTOR_ELT(R_points,i));
					  SET_VECTOR_ELT(R_ans, 2, R_I);
					  SET_VECTOR_ELT(R_ans, 3, R_S);
					  SET_VECTOR_ELT(R_ans, 4, R_jac);
					  SET_VECTOR_ELT(R_ans, 5, R_Iobs);
					  setVmatAttrib(&qlm, R_VmatNames, R_ans);

					  SET_VECTOR_ELT(R_ret, i, R_ans);
					  UNPROTECT(5);
			  }
          }
          UNPROTECT(4);
          return R_ret;
      }
}


/** \brief Quasi-Fisher score statistic
 *
 * @param x point
 * @param opt data struct
 * @param f value
 *
 */


double
ql_model_s::qfScoreStat(double *x, double *jac, double *score, double *qimat) {
	/* kriging statistics */
	glkm->intern_kriging(x);
	//printVector("krig.mean",glkm->krigr[0]->mean,&dx);
	//printVector("krig.var",glkm->krigr[0]->sig2,&dx);

	glkm->intern_jacobian(x,jac);
	//printMatrix("jac",jac,&dx,&nCov);

	intern_cvError(x);
	varMatrix(x,glkm->krigr[0]->sig2,qld->vmat);
	//printMatrix("vmat",qld->vmat,&nCov,&nCov);

	intern_quasiScore(jac,score);
	//printVector("score",score,&dx);

	intern_quasiInfo(jac,qimat);
	//printMatrix("Imat",qimat,&dx,&dx);

	return qfValue(score,qimat);
}


double
ql_model_s::intern_qfTrace(double *x) {
	int i=0;

	/* kriging */
	glkm->intern_kriging(x);
	glkm->intern_jacobian(x,jac);

	/* quasi-score */
	intern_cvError(x);
	varMatrix(x,glkm->krigr[0]->sig2,qld->vmat);
	intern_quasiScore(jac,score);

	/*! unfortunately transpose jac  and calculate matrix A */
	mat_trans(qld->Atmp,nCov,jac,dx,dx,nCov);
	solve_DSPTRS(qld->vmat,nCov,qld->Atmp,dx,&info);
	if(info > 0){
	  XERR(info,"Computation of variance of quasi-score failed.");
	  return R_NaN;
	}
	/* this actually computes the (average) trace of the variance of quasi-score */
	intern_varScore(qimat);
	double sum=0;
	for(; i < dx; ++i)
	  sum += qimat[i*dx+i];
	return sum/(double)dx;
}

double
ql_model_s::intern_qfVarStat(double *x) {
	/* kriging */
	glkm->intern_kriging(x);
	glkm->intern_jacobian(x,jac);

	/* quasi-score */
	intern_cvError(x);
	varMatrix(x,glkm->krigr[0]->sig2,qld->vmat);
	intern_quasiScore(jac,score);

	/*! unfortunately transpose jac  and calculate matrix A */
	mat_trans(qld->Atmp,nCov,jac,dx,dx,nCov);
	solve_DSPTRS(qld->vmat,nCov,qld->Atmp,dx,&info);
	if(info != 0){
	  XERR(info,"Computation of quasi-score variance failed.");
	  return R_NaN;
	}
	/* this actually computes the variance of quasi-score */
	intern_varScore(qimat);

	return qfValue(score,qimat);
}

/**
 * Modified quasi-deviance with
 * varS as weighting matrix instead of QI
 */
double
ql_model_s::qfValue(double *score, double *varS) {
	int i=0;
	MEMCPY(qld->score_work,score,dx);
	solve_DSPTRS(varS,dx,qld->score_work,ONE_ELEMENT,&info);
	if(info != 0){
	  XERR(info,"Computation of quasi-score statistic failed.");
	  return R_NaN;
	}

	double sum=0;
	for(; i < dx; ++i)
	  sum += qld->score_work[i]*score[i];
	return sum;
}


void
ql_model_s::intern_quasiInfo(double *jac, double *qimat) {
	/*! unfortunately transpose jac */
	mat_trans(qld->Atmp,nCov,jac,dx,dx,nCov);
	solve_DSPTRS(qld->vmat,nCov,qld->Atmp,dx,&info);
	if(info != 0){
	  XERR(info,"Computation of quasi-information matrix failed.");
	  return;
	}
	matmult(jac,dx,nCov,qld->Atmp,nCov,dx,qimat,&info);
	if(info)
	 WRR("`NaN` detected in `matmult`.")
}


void
ql_model_s::quasiScore(double *mean, double *jac, double *vmat, double *score) {
	int info=0, i=0;
	for(;i<nCov;i++)
	 qld->qtheta[i]=qld->obs[i]-mean[i];

	solve_DSPTRS(vmat,nCov,qld->qtheta,ONE_ELEMENT,&info);
	if(info != 0)
	  XWRR(info,"Computation of quasi-score (solving) failed: ");

	matmult(jac,dx,nCov,qld->qtheta,nCov,ONE_ELEMENT,score,&info);
	if(info)
	  WRR("`NaN` detected in matrix multiplication.")
}


/**  \brief Wrapper function: Score vector
 *
 * @param x     Parameter vector
 * @param data  Pointer data
 * @param score Score vector allocated
 */
void
wrap_intern_quasiScore(double *x, void *data, double *score) {
   ql_model qlm = (ql_model) data;
   GLKM glkm = qlm->glkm;
   ql_data qld = qlm->qld;

   krig_result krig_tmp = glkm->krigr[1];
   glkm->kriging(x,krig_tmp->mean,krig_tmp->sig2,krig_tmp->w);
   glkm->jacobian(x,krig_tmp->mean,qld->jactmp);

   if(!qld->qlopts.useSigma) {
     if(glkm->krigType && qld->qlopts.useCV)
       qlm->cvmod->cvError(x,glkm->km,krig_tmp->sig2);
     qlm->varMatrix(x,krig_tmp->sig2,qld->vmat_work);
   }

   qlm->quasiScore(krig_tmp->mean,qld->jactmp,qld->vmat_work,score);

}

/** \brief Wrapper function: intern_kriging
 *
 *     Comment: Kriging result is constructed and destroyed
 *
 * @param x point
 * @param data data pointer
 * @param mean kriging mean vector
 */
void
wrap_intern_kriging(double *x, void *data, double *mean) {
  GLKM glkm = (GLKM) data;
  krig_result krig_tmp = glkm->krigr[1];
  glkm->kriging(x,mean,krig_tmp->sig2,krig_tmp->w);
}

void
ql_model_s::varMatrix(double *x, double *s, double *vmat) {
	/* kriging prediction variances */
	if(qld->qlopts.varType == KRIG) {
	   if(varkm == NULL)
		 XERR(1,"Null pointer to Kriging models for covariance interpolation.");
	   varkm->intern_kriging(x);
	   /*! Merge to matrix*/
	   chol2var(varkm->krigr[0]->mean,vmat,nCov,qld->workx);
	   //printVector("varkm->mean",varkm->krigr[0]->mean,&nCov);
	   //printMatrix("vmat",vmat,&nCov,&nCov);

	   if(glkm->krigType) {
		   //printVector("s",s,&nCov);
		   add2diag(vmat,nCov,s);
	   }
	} else if(glkm->krigType) {
		 addVar(s,nCov,vmat,qld->work);
	}
}

void
ql_model_s::intern_quasiObs(double *x, double *score, double *qiobs) {
   fdJac(x,dx,score,dx,qiobs,&wrap_intern_quasiScore,(void*) this,FD_EPS,ONE_ELEMENT,&info);
   if(info)
	 WRR("`NaN` values detected in `fdJac`.")
}


/** @brief Allocate storage for a kriging model
 *
 * @param Xmat
 * @param lx
 * @param dx
 * @param trend
 * @param type
 * @return
 */

void
krig_model_s::alloc() {
  trend = cov.trend;
  switch (trend) {
    case 0: fddim = 1;
            break;
    case 1: fddim = dx+1;
            break;
    case 2: fddim = (dx+1)*(dx+2)/2;
            break;
    default: error(_("Invalid drift term specified.\n"));
             break;
  }
  CALLOCX(s0,lx,double);
  CALLOCX(f0,fddim,double);
  CALLOCX(Fmat,lx * fddim, double);
  CALLOCX(Cmat,lx * lx, double);

  /** allocate for inverse covariance matrix */
  if(krigtype == VARIANCE) {
      CALLOCX(Cinv,lx * lx, double);
      CALLOCX(X1mat,lx * fddim, double);
      CALLOCX(Qmat,fddim * fddim, double);
      CALLOCX(dw,lx + fddim, double);
  } else {
      CALLOCX(dw,lx + fddim, double);
  }
}

void
krig_model_s::setup(double *_data)
{
	  alloc();

	  if(cpy){
	   	CALLOCX(data,lx,double);
	   	MEMCPY(data,_data,lx);
	  } else {
	    data = _data;
	  }

	  // trend matrix
	  Fmatrix(Xmat,Fmat,lx,dx,trend);
	  // REML covariance matrix
	  if( intern_covMatrix(Xmat,dx,lx,Cmat,&cov) )
	    ERR("NaN` detected in covariance matrix computation.")

	  // init matrices
	  int info=0;
	  if(krigtype == VARIANCE) {
		 MEMCPY(Cinv,Cmat,lx*lx);
		 invMatrix(Cinv,lx,&info);
		 if(info)
		   XERR(info,"Setting up kriging model failed due to inversion error of covariance matrix.")
		 /* store matrices for solving kriging equations */
		 matmult(Cinv,lx,lx,Fmat,lx,fddim, X1mat,&info);
		 matmult_trans(Fmat,lx,fddim,X1mat,lx,fddim, Qmat,&info);
		 if(info){
			 WRR("`NaN` detected in matrix multiplication.")
		 }
		 /* init dual kriging equations */
		 intern_dualKrigingWeights(Cmat,lx,Fmat,fddim,data,dw);
	 } else {
		 /* init dual kriging equations */
		 intern_dualKrigingWeights(Cmat,lx,Fmat,fddim,data,dw);
	 }
}

void
cv_model_s::set(SEXP R_Xmat, SEXP R_data, SEXP R_cm) {
	int i=0, j=0, k=0, l=0, m=0,
		lx=0, len=0;

	int *id=0,
	    *dim=GET_DIMS(R_Xmat);

	double tmp=0;
	double *Xmatp=NULL,
		  **datap=NULL;

	np  = dim[0];
	dx  = dim[1];

	nc  = length(R_cm);
	fnc = 1.0/(nc*(nc-1));
	Xmat = REAL(R_Xmat);
	nCov = LENGTH(VECTOR_ELT(R_cm,0));
    errType = !std::strcmp("max",translateChar(asChar(getAttrib(R_cm,install("type")))));

	/* container of CV models*/
	if( (cm = new(std::nothrow) GLKM[nc]) == NULL)
		MEM_ERR(1,GLKM);

	CALLOCX(s2,nCov, double);
	CALLOCX(ybar,nCov, double);
	CALLOCX(ytil,nc*nCov, double);
	CALLOCX(y0,nCov,double);

	SEXP R_covList, R_id;
	/* construct new CV model data */
	for(k = 0; k < nc; k++) {
		R_covList = VECTOR_ELT(R_cm,k);
		R_id = AS_INTEGER(getAttrib(R_covList,install("id")));
		len = LENGTH(R_id);
		id = INTEGER(R_id);
		lx = np - len;

		CALLOCX(datap,nCov,double*);
		CALLOCX(Xmatp,lx*dx,double);

		for(j = 0; j < nCov; j++)
		 CALLOCX(datap[j],lx,double);

		// points
		for(m = 0, i = 0; i < np; i++) {
			 Rboolean isin = FALSE;
			 for(l = 0; l < len; l++) {
				if(i == (id[l]-1)) {
				  isin = TRUE;
				  break;
				}
			 }
			 if(isin) continue;
			 for(j = 0; j < nCov; j++) {
			   tmp = REAL(AS_NUMERIC(VECTOR_ELT(AS_LIST(R_data),j+dx)))[i];
			   if (ISNAN(tmp) || ISNA(tmp) || !R_finite(tmp) )
				   WRR("`NaN` detected in data vector.")
			   datap[j][m] = tmp;
			 }
			 for(j = 0; j < dx; j++)
			   Xmatp[lx*j+m] = Xmat[np*j+i];
			 ++m;
		}
		if( (cm[k] = new(std::nothrow) glkrig_models_s(R_covList,Xmatp,datap,lx,dx,krigType))==NULL)
			MEM_ERR(1,glkrig_models_s);

		// free memory
		for(j = 0; j < nCov; j++)
	      FREE(datap[j]);

		FREE(datap);
		FREE(Xmatp);
	}
}

void
cv_model_s::cvError(double *x, krig_model *km, double *cv) {
	int i = 0, k = 0;
	double yhat = 0;

	for(k = 0; k < nCov; ++k) {
		s2[k] = ybar[k] = 0;
		km[k]->dualKriging(x,ONE_ELEMENT,y0+k);
	}
	GLKM cmp=0;
	krig_model *kmp=0;

	for(i = 0; i < nc; ++i) {
			cmp = cm[i];
			kmp = cmp->km;
			for(k = 0; k < nCov; ++k) {
				kmp[k]->dualKriging(x,ONE_ELEMENT,&yhat);
				ytil[k*nc+i] = nc*y0[k] - (nc-1)*yhat;
				ybar[k] = ybar[k] + ytil[k*nc+i];
				if (!R_finite(ybar[k]) || ISNA(ybar[k]) || ISNAN(ybar[k])){
				  WRR("`NaN` detected in cross-validation errors.");
			      break;
				}
			}
	}
	for(i = 0; i < nc; ++i) {
		for(k = 0; k < nCov; ++k)
			s2[k] = s2[k] + SQR(ytil[k*nc+i]-ybar[k]/nc);
	}
	if(errType) {
	 for(k = 0; k < nCov; ++k) {
	     tmp = fnc*s2[k];
		 cv[k] = MAX(tmp,cv[k]);
		 if (!R_finite(cv[k]) || ISNA(cv[k]) || ISNAN(cv[k])) {
			 WRR("`NaN` detected in cross-validation errors.");
		     break;
		 }
	 }
	} else {
		for(k = 0; k < nCov; ++k){
		 cv[k] = fnc*s2[k];
		 if (!R_finite(cv[k]) || ISNA(cv[k]) || ISNAN(cv[k])) {
		   WRR("`NaN` detected in cross-validation errors.");
		   break;
		 }
		}
	}

}


inline void
glkrig_models::jacobian(double *x, double *mean, double *jac) {
	fdJac(x,dx,mean,nCov,jac,&wrap_intern_kriging,(void*)this,FD_EPS,ZERO_ELEMENT, &info);
	if(info)
	 WRR("`NaN` detected in `fdJac`.")
}

/**
 * Kriging a single point,
 * for all nCov statistics
 */
inline void
glkrig_models::kriging(double *x, double *m, double *s, double *w) {
	if(krigType) {
		for(int k=0;k < nCov; k++)
			km[k]->univarKriging(x,ONE_ELEMENT,m+k,s+k,w+k);
	} else {
		for(int k=0;k < nCov; k++)
			km[k]->dualKriging(x,ONE_ELEMENT,m+k);
	}
}


/**
 * Dual kriging a single point, single statistic
 */
void
krig_model_s::dualKriging(double *x, int nx, double *m) {
	intern_covVector(Xmat,dx,lx,x,nx,s0,&cov);
	trendfunc(x,nx,dx,f0,cov.trend);
	/** weights do not change */
	intern_dualKrigingPrediction(s0,f0,fddim,dw,lx,m);
}

/** \brief Internal function: Kriging
 * Chilès, J.-P. and Delfiner, P. (2008) Frontmatter, in Geostatistics:
 *              Modeling Spatial Uncertainty, John Wiley & Sons, Inc., Hoboken, NJ, USA.
 *              page: 169
 *
 * @param Xmat sample matrix
 * @param data observed data vector
 * @param Fmat trend matrix
 * @param fddim columns of trend matrix
 * @param Cinv inverse covariance matrix (stored)
 * @param point prediction point x0 as matrix
 * @param Npoints rows of point
 * @param mean predicted (kriging) mean
 * @param sigma2 prediction (kriging) variance
 * @param w kriging weights
 * @param cov covariance model
 */
void
krig_model_s::univarKriging(double *x, int nx, double *m, double *s, double *w) {
	 /*!   X1 = C^{-1} F
      *     Q = F^{t} X1
	  */

	/*! calculates the covariance vector: sigma0 */
	intern_covVector(Xmat,dx,lx,x,nx,s0,&cov);
	//printVector("inter_covVector",sigma0,lx);

	/*! calculates the trend vector: f0 */
	trendfunc(x,nx,dx,f0,cov.trend);
	//printVector("inter_trendfunc",f0,fddim);

	*s = cov.cf(&cov,&ZERO_ELEMENT);
	//printMatrix("Qmat",Qmat,fddim,fddim);

	/*! calculate kriging weights and kriging variance */
	solveKriging(Cinv,Fmat,X1mat,Qmat,lx,fddim,s,w,s0,f0);

	/* manipulate Kriging MSE */
	//*s += cov.nugget;

	/*! kriging mean */
	*m = 0.0;
	for (int j=0; j<lx; j++)
	  *m += w[j] * data[j];
}



#endif /* KRIGING_CC_ */
