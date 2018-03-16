/**
 * @file        covariance.h
 * @date        03.07.2012
 * @author      M. Baaske
 * @brief       Covariance calculations
 *
 * Explanation:
 *
 */

#ifndef COVARIANCE_H_
#define COVARIANCE_H_

#include "error.h"
#include "auxil.h"

struct cov_model_s;

typedef double cov_param;
typedef double (*cov_func)(cov_model_s *, double *);

double covMatern(cov_model_s *cov, double *h);
double covSirfk(cov_model_s *cov, double *h);
double covPowExp2(cov_model_s *cov, double *h);
double covExp(cov_model_s *cov, double *h);

///-> no Calloc(,) here!
typedef struct cov_model_s {
    cov_func cf;
	cov_param *param;

    int npar,
	    trend,
		ln; 					// length of nuggfix

    double nugget,				// global nugget (estimated by reml)
	      *nuggfix;				// local nugget variances

    cov_model_s(SEXP R_Cov);

    ~cov_model_s() {
    	FREE(param);
    	FREE(nuggfix);
    }

} cov_model;

#ifdef  __cplusplus
extern "C" {
#endif

SEXP getListElement (SEXP list, const char *str);

void matern(double *scale, double *nu,double *rho, double *h);
void sirfk(double *h, int *dim, double *_alpha, double *_scale, double *eps, double *out);

SEXP covValue(SEXP R_h, SEXP R_Cov);
SEXP covFct(SEXP R_Xmat, SEXP R_Ymat, SEXP R_Cov);
SEXP covVector (SEXP R_Xmat, SEXP R_Yvec);
SEXP covMatrix( SEXP R_Xmat, SEXP R_cov);

SEXP Fmat(SEXP R_Xmat, SEXP R_trend);
SEXP Pmat(SEXP R_Fmat);

int intern_covMatrix(double *x, int dx, int lx, double *z, cov_model *cov);
int intern_covVector(double *x, int dx, int lx, double *y, int ly, double *z, cov_model* cov);

#ifdef  __cplusplus
}
#endif


#endif /* COVARIANCE_H_ */
