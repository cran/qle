/**
 * @file        covariance.cc
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       Covariance calculations
 *
 * Explanation:
 *
 */


#include "covariance.h"

#include <Rmath.h>

#define ZERO_ELEMENT 0
#define MIN_DISTANCE 1e-16
#define BESSELK_TOL 1e-16
#define MATERN_NU_TOL 100
#define MATERN_RESULT_TOL 1e-16

SEXP getListElement (SEXP list, const char *str) {
     SEXP elmt = R_NilValue;
     SEXP names = getAttrib(list, R_NamesSymbol);

     for (int i = 0; i < length(list); i++)
         if(std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
             elmt = VECTOR_ELT(list, i);
             break;
         }
     return elmt;
}

cov_model_s::cov_model_s(SEXP R_Cov) :
  	param(0),npar(0),
	trend(0),ln(0),nugget(0),nuggfix(0)
{
  	SEXP R_nuggfix, R_param;
  	int cvtype = asInteger(getListElement(R_Cov, "model"));

  	switch(cvtype) {
	   case 1:
		 cf = covSirfk;
		 break;
	   case 2:
		 cf = covMatern;
		 break;
	   case 3:
		 cf = covPowExp2;
		 break;
	   case 4:
	   	cf = covExp;
	   	break;
	   default: ERR("unknown covariance model type."); break;
	}
  	/* trend order */
  	trend = asInteger(getListElement(R_Cov, "trend"));

	/* model parameter */
	R_param = getListElement(R_Cov, "param");
	if(isNull(R_param) || !isNumeric(R_param))
	  ERR("Covariance parameters are `Null` or include non-numeric values.");

	npar = LENGTH(R_param);
	CALLOCX(param,npar,double);
	MEMCPY(param,REAL(R_param),npar);
	nugget = param[npar-1];

	/* local nugget variance (for each sampled location) */
	R_nuggfix = getListElement(R_Cov, "fix.nugget");

	if(!isNull(R_nuggfix)){
		if(!isNumeric(R_nuggfix))
		  ERR("The vector of local nugget variances should be numeric.");
		ln = LENGTH(R_nuggfix);
	    if(ln > 0){
			CALLOCX(nuggfix,ln,double);
			MEMCPY(nuggfix,REAL(R_nuggfix),ln);
		} else WRR("local nugget variances are replicated.")
	}

}


/*
  m = length(R_cp);
  printVector("param", cov->param, &m);
  m = length(R_cfp);
  printVector("nuggfix",cov->nuggfix,&cov->npoints);
  DEBUG_DUMP_VAR(cov->nugget,"%f");
  DEBUG_PRINT_VECTOR_R(R_nugg,"nugget");
*/


SEXP covVector (SEXP R_Xmat, SEXP R_Yvec, SEXP R_CovStruct) {
    int info=0,
    	*dimX = GET_DIMS(R_Xmat);

    int lx = dimX[0],
    	dx = dimX[1];

    if(!isVector(R_CovStruct))
	 ERR("Expected list with covariance parameters.\n");

    cov_model cov(R_CovStruct);

    SEXP R_v=R_NilValue;
    PROTECT(R_v = allocVector(REALSXP, lx));

    if( (info = intern_covVector(REAL(R_Xmat), dx, lx, REAL(R_Yvec), 1, REAL(R_v), &cov)) != NO_ERROR )
       LOG_WARNING(info," in `inter_covVector`.");

    UNPROTECT(1);
    return R_v;
}

int intern_covVector(double *x, int dx, int lx, double *y, int ly, double *z, cov_model* cov) {
	int have_na=0;
	double h=0, *px=x;
	cov_func cf=cov->cf;

	for (int i = 0; i < lx; i++, px++) {
		h = norm_2(px, y, lx, ly, dx);
		/* filter out nugget-effect component, that is, estimate
		 * the simulation variance free value of the statistics */
		z[i] = (*cf)(cov,&h);
		if(!R_FINITE(z[i]))
		  { have_na=1; break; }

	    /* not filtering out the nugget-effect component, that is, kriging as an exact interpolator */
		//if( h < MIN_DISTANCE) {
		//   z[i] = (*cf)(cov,&zero) + cov->nugget + cov->nuggfix[i];
		//} else {
		//   z[i] = (*cf)(cov,&h);
		//}
	}
	return have_na;
}


/**
 * @brief
 *
 * @param R_X          Distance matrix
 * @param R_pars       Covariance parameters
 * @param R_nugg       Global nugget value
 * @param R_nuggfix    Local nugget vector
 * @return
 */
SEXP covMatrix(SEXP R_Xmat, SEXP R_cov ) {
  int info = 0,
	  *dimX = GET_DIMS(R_Xmat);
  int lx = dimX[0],
	  dx = dimX[1];

  cov_model cov(R_cov);
  SEXP R_C = R_NilValue;
  PROTECT(R_C = allocMatrix(REALSXP,lx,lx));

  if( (info = intern_covMatrix(REAL(R_Xmat),dx,lx,REAL(R_C),&cov)) != NO_ERROR )
    LOG_WARNING(info, " in `intern_covMatrix`.");

  UNPROTECT(1);
  return R_C;
}


int intern_norm(double *x, int lx, int dx, double *z) {
	int have_na =0;
    double h=0;
    for (int i = 0; i < lx; i++) {
        for (int j = 0; j < i; j++) {
            h = norm2(x, lx, x, lx, dx, i, j);
            if(!R_FINITE(h))
             { have_na=1; break; }
            z[j+lx*i] = z[i+lx*j] = h;
        }
        z[i+lx*i] = 0;
    }
    return have_na;
}

int intern_covMatrix(double *x, int dx, int lx, double *z, cov_model *cov) {
    int have_na =0;
    double val = 0,
             h = 0,
           zero= 0;

    cov_func cf = cov->cf;
    if(cov->nuggfix) {
    	  if(cov->ln != lx){
    		  Rprintf("length local nugget variances: %d, locations: %d \n",cov->ln,lx);
    		  ERR("Length of local nugget variances does not match number of locations.")
    	  }
		  for (int i = 0; i < lx; i++) {
			for (int j = 0; j < i; j++) {
				h = norm2(x, lx, x, lx, dx, i, j);
				val = (*cf)(cov,&h);
				if (!R_FINITE(val) )
				  { have_na=1; break; }
				z[j + lx * i] = z[i + lx * j] = val;
			}
			z[i + lx *i] =(*cf)(cov,&zero) + cov->nugget + cov->nuggfix[i]; //diagonal terms
		  }
    } else {
    	for (int i = 0; i < lx; i++) {
			for (int j = 0; j < i; j++) {
				h = norm2(x, lx, x, lx, dx, i, j);
				val = (*cf)(cov,&h);
				if(!R_FINITE(val))
				 { have_na=1; break; }
				z[j + lx * i] = z[i + lx * j] = val;
			}
			z[i + lx *i] =(*cf)(cov,&zero) + cov->nugget;
		}
    }
    //printMatrix("covMatrix: ",z,&lx,&lx);
    return have_na;
}

SEXP Pmat(SEXP R_Fmat){
    int info = 0;
    int *dimX = GET_DIMS(R_Fmat);
    int lx = dimX[0];
    int fddim = dimX[1];

    SEXP R_Pmat = R_NilValue;
    PROTECT(R_Pmat = allocMatrix(REALSXP, lx, lx-fddim));   	// nullspace matrix

    qr_data qr = qr_dgeqrf(REAL(R_Fmat), lx, fddim, info);
    nullSpaceMat(qr,REAL(R_Pmat), info);
    qrFree(qr);

    UNPROTECT(1);
    return R_Pmat;
}

SEXP Fmat(SEXP R_Xmat, SEXP R_trend) {
    if(!isMatrix(R_Xmat))
     ERR("Sample data must be a matrix order.");

    int *dimX = GET_DIMS(R_Xmat);
    int lx = dimX[0],
    	dx = dimX[1];

    /* Size F matrix: stored columns first
	*
	* constant trend:  size(x) x 1
	* linear trend  :  size(x) x (dim + 1)
	* quadratic trend: size(x) x (dim+1)*(dim+2)/2
	*
    */

    int fddim = 0;
    int trend = asInteger(R_trend);
    switch (trend) {
		case 0: fddim = 1; break;
		case 1: fddim = dx+1; break;
		case 2: fddim = (dx+1)*(dx+2)/2; break;
		default:
		  error(_("Invalid drift term specified."));
		  break;
    }

    SEXP R_Fmat = R_NilValue;
    PROTECT(R_Fmat=allocMatrix(REALSXP, lx, fddim));

    /* construct Trend matrix */
    Fmatrix(REAL(R_Xmat),REAL(R_Fmat),lx,dx,trend);

    UNPROTECT(1);
    return R_Fmat;
}


SEXP covFct(SEXP R_Xmat, SEXP R_Ymat, SEXP R_covStruct) {
    if(!isMatrix(R_Xmat) || !isMatrix(R_Ymat))
       error(_("Expected matrix as first input argument!"));
    int *dimX = GET_DIMS(R_Xmat);
    int lx = dimX[0],
        dx = dimX[1];

    if(!isVector(R_covStruct))
       ERR("Expected covariance model as list argument in `covFct`");

    if( lx != GET_DIMS(R_Ymat)[0] || dx != GET_DIMS(R_Ymat)[1]) {
       ERR("Matrix dimensions do not match in `covFct`.");
    }
    cov_model cov(R_covStruct);

    SEXP R_v = R_NilValue;
    PROTECT(R_v = allocVector(REALSXP, lx));

    double *cv = REAL(R_v);
    double *x = REAL(R_Xmat);
    double *y = REAL(R_Ymat);

    int lxdx = lx *dx;
    double *dhx = (double*) R_alloc( (size_t) lxdx, sizeof(double));
    dist_X1_X2(x,lx,dx,y,lx,dhx);

    double h = 0;
    for(int k=0; k < lx; k++, cv++) {
      h = norm_x_2(dhx,lx,dx,k);
      *cv = (cov.cf)(&cov,&h);
    }

    UNPROTECT(1);
    return R_v;
}

SEXP covValue(SEXP R_h, SEXP R_covStruct ) {
    if(!isVector(R_covStruct))
       ERR("Expected argument of type list in `covValue`.");

    cov_model cov(R_covStruct);
    return ScalarReal((cov.cf)(&cov,REAL(AS_NUMERIC(R_h))));
}

void sirfk(double *_h, int *_dim, double *_alpha, double *_scale, double *_eps, double *out) {
  double alpha = *_alpha,
         h = *_h,
         eps = *_eps,
         scale = *_scale,
         a = alpha-std::floor(alpha),
         b = std::ceil(alpha)-alpha,
         delta =( a < b) ? a : b;

  int d = *_dim, fac = 1;

  if(h < MIN_DISTANCE) {
      *out = ZERO_ELEMENT;
  } else {
   //alpha integer
   if( delta < eps ) {
      alpha = ( a < b) ? std::floor(alpha) : std::ceil(alpha);

      for(int i=1;i<=alpha;i++) fac=fac*i;
      *out = scale * std::pow(2,-2*alpha+1)*std::pow(-1,alpha-1)*pow(M_PI,d/2)/(gammafn( (double) alpha+d/2)*fac)*std::pow(h,2*alpha)*std::log(h);
   }
   //alpha not an integer
   else {
      *out = scale * std::pow(2,-2*alpha)*std::pow(M_PI,d/2)*gammafn(-alpha)/gammafn( (double) alpha+d/2) * std::pow(h,2*alpha);
   }
  }
}

double covSirfk(cov_model *cov, double *hv){
	if(*hv < MIN_DISTANCE)
	   return ZERO_ELEMENT;

	double h = *hv,
		   alpha = cov->param[1],
           a = alpha-std::floor(alpha),
           b = std::ceil(alpha)-alpha;

    if( MIN(a,b) < DBL_EPSILON ) {
         alpha = ( a < b) ? std::floor(alpha) : std::ceil(alpha);                // force alpha as integer value
         return cov->param[0] * std::pow(-1,alpha+1)*std::pow(h,2*alpha)*std::log(h);
    }
    return cov->param[0] * std::pow(-1,floor(alpha)+1)*std::pow(h,2*alpha); 	// alpha not integer
}

void matern(double *scale, double *nu,double *rho, double *h) {
    double nuTol = *nu < MATERN_NU_TOL ? *nu : MATERN_NU_TOL,
    y = *h * std::sqrt(2.0 * nuTol)/ *rho;

    static double nuOld=R_PosInf;
    static double lgamma;

    if(y <= BESSELK_TOL)
        *h = 1.0;
    else {
        if (*nu != nuOld) {
            nuOld = *nu;
            lgamma = lgammafn(*nu);
        }
       *h = *scale * (2.0* std::exp(*nu * std::log(0.5 * y) - lgamma + std::log(bessel_k(y, *nu, 2.0)) - y));
    }
}


double covMatern(cov_model *cov, double *h) {
    double nu = cov->param[1],
          rho = cov->param[2];

    double nuTol = nu < MATERN_NU_TOL ? nu : MATERN_NU_TOL,
    y = *h * std::sqrt(2.0 * nuTol)/rho;
    if(!R_FINITE(y)){
    	WRR("`NaN` detected in `covMatern`.")
    	return R_NaN;
    }

    static double nuOld=R_PosInf;
    static double lgamma;

    double ret;
    if(y <= BESSELK_TOL)
      ret = 1.0;
    else {
      if (nu != nuOld) {
        nuOld = nu;
        lgamma = lgammafn(nu);
      }
     ret = 2.0 * std::exp(nu * log(0.5 * y) - lgamma + std::log(bessel_k(y, nu, 2.0)) - y);
    }
    return cov->param[0] * ret;
}

double covPowExp2(cov_model *cov, double *h) {
    if(*h<MIN_DISTANCE)
      return cov->param[0];

    double phi = cov->param[1];
    double kappa = cov->param[2];

    return( cov->param[0] * std::exp(-std::pow(*h/phi,kappa)) );
}

double covExp(cov_model *cov, double *h) {
	 if(*h<MIN_DISTANCE)
	   return cov->param[0];

	 return (cov->param[0] * std::exp( -( *h/cov->param[1] ) ) );
}
