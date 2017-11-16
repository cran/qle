/**
 * @file        kriging.h
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions for kriging
 *
 *
 */

#ifndef KRIGING_H_
#define KRIGING_H_

#include "error.h"
#include "covariance.h"

/** typedefs */
typedef enum { MEAN = 0, KRIG = 1} var_type;
typedef enum { UNKNOWN = -1, DUAL = 0, VARIANCE = 1, BOTH = 2} krig_type;

typedef struct krig_result_s {
  double *mean,
         *sig2,
            *w;

  krig_type krigType;

  krig_result_s() :
  	  mean(0), sig2(0), w(0), krigType(VARIANCE)
  {}

  krig_result_s(int nCov, int lx, krig_type type) :
	  mean(0), sig2(0), w(0), krigType(type)
  {
	alloc(nCov,lx);
  }

  ~krig_result_s() {
	 FREE(w)
	 FREE(mean)
	 FREE(sig2)
  }

  void alloc(int nCov, int lx) {
	  if(nCov==0 || lx==0)
	  	 error(_("number of covariance models or sample points must not be zero!"));

	  CALLOCX(w,nCov*lx,double);     // kriging weights as matrix: nCov rows, lx columns
	  CALLOCX(mean,nCov,double);     // kriging means (at single point)
	  if(krigType)
	     CALLOCX(sig2,nCov,double);  // kriging variances (at single point) and nCov covariances
  }

} krig_result_t, *krig_result;

typedef struct ql_options_s {
  var_type varType;
  Rboolean useCV, useSigma;

  ql_options_s(SEXP R_opts)  {
	  useCV    = (Rboolean) asLogical(getListElement(R_opts, "useCV"));
	  useSigma = (Rboolean) asLogical(getListElement(R_opts, "useSigma"));

	  SEXP R_varType = getListElement(R_opts, "varType");
	  if(!isString(R_varType))
	    error(_("Variance interpolation type should be a string: varType"));
	  const char *var_type = translateChar(asChar(R_varType));
	  if( !strcmp("kriging",var_type) )  {
	   varType = KRIG;
	  } else
	   varType = MEAN;
  }
} ql_options_t, *ql_options;


typedef struct krig_model_s {
  int lx,      /* amount of sampled points */
	  dx;      /* dimension */

  Rboolean cpy;	 /* if 'data' should be  copied */
  krig_type krigtype;

  double *Xmat,	 /* design (always only pointer) */
  	  	 *data,  /* stats vector */
         *Fmat,  /* Trend matrix */
         *Cmat,  /* Covariance matrix */
         *Cinv,  /* Inverse Covariance */
         *X1mat, /* X1 matrix, solving kriging equations */
         *Qmat,  /* Q matri, solving kriging equations */
         *dw;    /* Dual kriging weights */

  double *f0,*s0; /* trend and covariance vector work */
  int   trend,   /* trend order */
        fddim;   /* columns of trend matrix F*/

  cov_model cov;

  krig_model_s(SEXP _R_cov, double *_Xmat, double *_data,
		  int _lx, int _dx, krig_type _type,  Rboolean _cpy = TRUE 	/* copy _data */ ) :
	 lx(_lx), dx(_dx), cpy(_cpy), krigtype(_type), Xmat(_Xmat),
	 data(0),Fmat(0),Cmat(0),Cinv(0),X1mat(0),Qmat(0),dw(0),f0(0),s0(0),trend(0),fddim(0),
	 cov(_R_cov)
  {
	 setup(_data);
  }

  ~krig_model_s() {
	  if(krigtype) {
	      FREE(Cinv)
	      FREE(X1mat)
	      FREE(Qmat)
	      FREE(dw)
	  } else {
	      FREE(dw)
	  }
	  FREE(s0)
	  FREE(f0)
	  FREE(Fmat)
  	  FREE(Cmat)
  	  if(cpy)
  		FREE(data)
  }

  void alloc();

  void setup(double *data);

  void dualKriging(double *x, int nx, double *m);

  void univarKriging(double *x, int nx, double *m, double *s, double *w);


} krig_model_t, *krig_model;


typedef struct glkrig_models_s
{
  double *Xmat;

  int lx,         /* rows of sample matrix */
      dx,         /* columns, dimension of sample points */
      nCov;

  Rboolean cpy;
  krig_type krigType;
  int info;

  krig_model *km;
  krig_result krigr[2];   /* working storage kriging: 0: default, 1: temporary */

  glkrig_models_s(SEXP _R_covList, SEXP _R_Xmat, SEXP _R_data, SEXP R_krigType, int idx = 0, Rboolean _cpy = TRUE) :
	  Xmat(0), lx(0), dx(0), nCov(0), cpy(_cpy), info(0)
  {
	  int *dim = GET_DIMS(_R_Xmat);
	  lx = dim[0];
	  dx = dim[1];

	  if(isNull(_R_covList))
		ERR("Covariance model is not set (NULL).");
	  nCov = LENGTH(_R_covList);

	  if(isNull(R_krigType))
		 krigType = VARIANCE;
	  else {
		  if(!isString(R_krigType))
	  		error(_("Kriging type estimation should be a string"), "R_krigType");
		  const char *krig_type = translateChar(asChar(R_krigType));
		  if( !strcmp("dual",krig_type) )
			  krigType = DUAL;
		  else if( !strcmp("var",krig_type) )
			  krigType = VARIANCE;
		  else
			  ERR("unknown kriging type: choose `var` or `dual`.");
	  }

	  /* kriging results storage */
	  for(int i=0; i<2; i++)
	  {
		if( (krigr[i] = new krig_result_s(nCov,lx,krigType)) == NULL)
		  MEM_ERR(1,krig_result_s)
	  }

	  if(cpy) {
		  CALLOCX(Xmat,lx*dx,double);
		  MEMCPY(Xmat,REAL(_R_Xmat),lx*dx);
	  } else {
		  Xmat = REAL(_R_Xmat);
	  }

	  if( (km =  new(std::nothrow) krig_model[nCov]) == NULL)
		  MEM_ERR(1,krig_model)

	  for(int i=0; i<nCov; i++)
	  {
		 if( ( km[i] = new(std::nothrow) krig_model_s(VECTOR_ELT(_R_covList,i),
				 	 	 	 	  Xmat, REAL(AS_NUMERIC(VECTOR_ELT(AS_LIST(_R_data),i+idx))),
								  lx,dx,krigType,cpy) ) == NULL )
			 MEM_ERR(1,krig_model_s)
	  }

  }

  glkrig_models_s(SEXP _R_covList, double *_Xmat, double **_data, int _lx, int _dx,
		  	  	  	 krig_type _krigType, int idx = 0, Rboolean _cpy = TRUE) :
		 Xmat(0), lx(_lx), dx(_dx), nCov(LENGTH(_R_covList)), cpy(_cpy), krigType(_krigType), info(0)
  {

  	  /* kriging results storage */
  	  for(int i=0; i<2; i++)
  	  {
  		if( (krigr[i] = new krig_result_s(nCov,lx,krigType)) == NULL)
  	  	  MEM_ERR(1,krig_result_s)
  	  }

  	  if(cpy) {
  		 CALLOCX(Xmat,lx*dx,double);
  		 MEMCPY(Xmat,_Xmat,lx*dx);
  	  } else {
  		Xmat = _Xmat;
  	  }

  	  if( (km =  new(std::nothrow) krig_model[nCov]) == NULL)
  		  MEM_ERR(1,krig_model)

  	  for(int i=0; i<nCov; i++)
  	    if( (km[i] = new(std::nothrow) krig_model_s(VECTOR_ELT(_R_covList,i),Xmat,
  	    		 	 	 _data[i+idx],lx,dx,krigType,cpy)) == NULL)
  	    	MEM_ERR(1,krig_model_s)
  }

  ~glkrig_models_s() {
	 int k=0;
	  if(km != NULL) {
	  for(k=0; k<nCov;k++)
	    DELETE(km[k])
	  DELETE_ARR(km)
	 }
	 for(k=0; k<2; k++)
	   DELETE(krigr[k])
	 if(cpy)
	   FREE(Xmat)
   }

  void kriging(double *x, double *m, double *s, double *w);
  inline void intern_kriging(double *x) {  kriging(x, krigr[0]->mean, krigr[0]->sig2, krigr[0]->w);}

  void jacobian(double *x, double *mean, double *jac);
  inline void intern_jacobian(double *x, double *jac) { jacobian(x, krigr[0]->mean, jac);  }

} glkrig_models, *GLKM;

typedef struct cv_model_s {
	int nc,			/* number of CV models, each made up of nCov statistics */
		np,dx,	    /* overall number and dimension of sample point(s) 	*/
		nCov,
		errType;
	double fnc, tmp;

	GLKM *cm;
	krig_type krigType;
	double *Xmat, *s2, *ybar, *ytil, *y0;

	cv_model_s(SEXP R_Xmat, SEXP R_data, SEXP R_cm) :
	 nc(0), np(0), dx(0), nCov(0), errType(0), fnc(0), tmp(0),
	 cm(0), krigType(VARIANCE),
	 Xmat(0), s2(0), ybar(0), ytil(0), y0(0)
	{
		set(R_Xmat, R_data, R_cm);
	}

	~cv_model_s() {
	  if(cm != NULL) {
		 for(int i=0; i<nc; i++)
		   DELETE(cm[i])
		 DELETE_ARR(cm)
	  }
	  FREE(y0)
	  FREE(s2)
	  FREE(ybar)
	  FREE(ytil)
	}

	void set(SEXP R_Xmat, SEXP R_data, SEXP R_cm);

	void cvError(double *x, krig_model *km, double *cv);

} cv_model_t, *CVKM;

/*! \todo:  kriging models for var(T(X)) */
typedef struct ql_data_s {

  // always to be allocated
  double *obs,          /* Statistics of reference parameter */
         *qtheta,       /* work: nCov */
		 *vmat,
		 *vmat_work,    /* variance matrix of statistics */
		 *score_work,	/* working score vector */
         *work,	        /* copy of diagonal terms of variance matrix of statistics */
		 *workx,		/* for kriging the variance matrix */
		 *jactmp, 		/* transpose of jac or simply temporary storage */
		 *Atmp, 	    /* temporary working matrix */
  	  	 *tmp;

  ql_options_t qlopts;

  int dx,lx,nCov; /* dimension parameter, number of samples, number of statistics */

  // alloc and init
  ql_data_s(SEXP R_obs, SEXP R_Vmat, SEXP R_qlopts, int _nCov, int *_dims) :
		  obs(0), qtheta(0), vmat(0), vmat_work(0), score_work(0), work(0),
		  workx(0), jactmp(0), Atmp(0), tmp(0), qlopts(R_qlopts),
		  dx(_dims[1]), lx(_dims[0]), nCov(_nCov)
  {
	 int nCov2 = SQR(nCov);

	 CALLOCX(qtheta,nCov,double);
	 CALLOCX(Atmp,dx*nCov,double);
	 CALLOCX(obs,nCov,double);
	 CALLOCX(tmp,nCov,double);
	 CALLOCX(jactmp,dx*nCov,double);
	 CALLOCX(vmat,nCov2,double);
	 CALLOCX(vmat_work,nCov2,double);
	 CALLOCX(score_work,dx,double);

	 for(int i=0;i<dx*nCov;i++)
		jactmp[i]=Atmp[i]=0;

	 double *_obs = REAL(AS_NUMERIC(R_obs));
	 for(int i=0;i<nCov;i++) {
		 tmp[i]=obs[i]=_obs[i];
	 }

	 if(!isNull(R_Vmat))
	 {
		if(!isMatrix(R_Vmat))
		  error(_("'%s' Variance is not of type matrix"), "R_Vmat");

		 if(qlopts.useSigma || qlopts.varType == MEAN) {
		    double *_vmat = REAL(R_Vmat);
		    MEMCPY(vmat,_vmat,nCov2);
		    MEMCPY(vmat_work,_vmat,nCov2);

		    if(qlopts.varType == MEAN) {
		      CALLOCX(work,nCov,double);
		      splitDiag(vmat_work,nCov,work);
		    }
		 }

	 } else {
		 CALLOCX(workx,nCov2,double);
	 }

  };

  ~ql_data_s() {
	  FREE(jactmp)
	  FREE(tmp)
	  FREE(obs)
	  FREE(Atmp)
	  FREE(qtheta)
	  FREE(vmat)
	  FREE(vmat_work)
	  FREE(score_work)
	  FREE(work)
	  FREE(workx)
  };

} ql_data_t, *ql_data;



typedef struct ql_model_s {
	//ql_data_t qld;

	GLKM glkm;  			/* Kriging models for statistics */
	GLKM varkm;  			/* Kriging models for variance matrices */
	CVKM cvmod;				/* CV models (Cross-Validation kriging models) */
	ql_data qld;

	double *qimat,  /* expected Quasi-Information */
		   *score,  /* quasi score vector */
		   *jac,    /* gradient statistics */
		   *qiobs;  /* observed Quasi-Information = dS/dtheta by FD */

	double *lower, *upper;   /* bounds on the parameter space */

	int dx,dxdx,lx,nCov,nCov2,info;

	// R_qsd, R_qlopts, R_Vmat, R_cm, (Rboolean) type
	ql_model_s(SEXP R_qsd, SEXP R_qlopts, SEXP R_Xmat, SEXP R_Vmat, SEXP R_cm, Rboolean _cpy) :
			glkm(NULL), varkm(NULL), cvmod(NULL), qld(NULL),
			qimat(0), score(0), jac(0), qiobs(0), lower(0), upper(0), info(0)
	{
		/* bounds */
		lower = REAL(getListElement( R_qsd, "lower" ));
		upper = REAL(getListElement( R_qsd, "upper" ));

		/* list of covariance models*/
		SEXP R_covT = getListElement( R_qsd, "covT" );
	    /* list of covariance models for variances*/
		SEXP R_covL = getListElement( R_qsd, "covL" );
		/* reference statistic values*/
		SEXP R_obs = getListElement( R_qsd, "obs" );
		/* ql data list */
		SEXP R_qldata = getListElement( R_qsd, "qldata" );
		SEXP R_krigType = getListElement( R_qsd, "krig.type" );

		nCov = LENGTH(R_covT);
		if( (qld = new(std::nothrow) ql_data_s(R_obs,R_Vmat,R_qlopts,nCov,GET_DIMS(R_Xmat))) == NULL)
		    MEM_ERR(1,ql_data_s)

		lx = qld->lx;
		dx = qld->dx;
		dxdx = SQR(dx);
		nCov2 = SQR(nCov);

		 if(!qld->qlopts.useSigma && qld->qlopts.varType)
		 {
			 if(!isNull(R_Vmat))
			   ERR("R_Vmat is not NULL but should be if using kriging approximation of variance.");
			 if(isNull(R_covL))
			   ERR("Covariance model for kriging variance matrix is not set (Null).");

			 int idx = dx+2*nCov;
			 if( (varkm = new(std::nothrow) glkrig_models_s(R_covL,R_Xmat,R_qldata,R_NilValue,idx,_cpy)) == NULL)
			   MEM_ERR(1,glkrig_models_s)
		 }
		 if(!isNull(R_cm) && qld->qlopts.useCV) {
			 if( (cvmod = new(std::nothrow) cv_model_s(R_Xmat,R_qldata,R_cm)) == NULL)
			  MEM_ERR(1,cv_model_s)
		 }

		 if( (glkm = new(std::nothrow) glkrig_models_s(R_covT,R_Xmat,R_qldata,R_krigType,dx,_cpy)) == NULL)
			MEM_ERR(1,glkrig_models_s)

		 CALLOCX(score,dx,double);
		 CALLOCX(qiobs,dxdx,double);
		 CALLOCX(qimat,dxdx,double);
		 CALLOCX(jac,dx*nCov,double);
	}

	~ql_model_s()
	{
	  DELETE(glkm)
	  DELETE(varkm)
	  DELETE(cvmod)
	  DELETE(qld)

	  FREE(qiobs)
	  FREE(qimat)
	  FREE(jac)
	  FREE(score)
	}

	inline void intern_cvError(double *x) {
		  if(glkm->krigType && qld->qlopts.useCV) {
		  	 cvmod->cvError(x,glkm->km,glkm->krigr[0]->sig2);
		  	 if(info)
		  	   WRR("`NaN` detected in cross-validation errors.")
		  }
    }

	  /* variance score vector, only for krigType `var` !!! */
	inline void intern_varScore(double *vars)
	{
   		 matmult_diag_sqrt(qld->Atmp,nCov,dx,glkm->krigr[0]->sig2);
		 matmult_trans(qld->Atmp,nCov,dx,qld->Atmp,nCov,dx,vars,&info);
		 if(info){
		   WRR("`NaN` detected in matrix multiplication.")
		 }

#if DEBUG
		printMatrix("Atmp (1)",qld->Atmp,&nCov,&dx);
		printVector("sig2",glkm->krigr[0]->sig2,&nCov);
		printMatrix("vars",vars,&dx,&dx);
#endif

    }

	inline void intern_quasiObs(double *x, double *score, double *qiobs);

	inline void intern_quasiScore(double *jac, double *score) {
		 quasiScore(glkm->krigr[0]->mean,jac,qld->vmat,score);
	}

	double intern_qfTrace(double *x);
    double intern_qfVarStat(double *x);

    inline double intern_qfScoreStat(double *x) {
	  	return qfScoreStat(x,jac,score,qimat);
	}

	void varMatrix(double *x, double *s, double *vmat);

	double qfValue(double *score, double *varS);
	double qfScoreStat(double *x, double *jac, double *score, double *qimat);

	void intern_quasiInfo(double *jac, double *qimat);
	void quasiScore(double *mean, double *jac, double *vmat, double *score);

	double intern_mahalValue(double *x);			//use constant (inverted) variance
	double intern_mahalVarTrace(double *x);			//use constant (inverted) variance
	double intern_mahalValue_theta(double *x);		// computing inverse variance


} ql_model_t, *ql_model;

//////////////////////////// wrapper ///////////////////////////////////////////////////////
void addVar(double *sig2, int nc, double *vmat, double *work);
void wrap_intern_kriging(double *x, void *data,  double *mean);
void wrap_intern_quasiScore(double *x, void *data, double *score);
void setVmatAttrib(ql_model qlm, SEXP R_VmatNames, SEXP R_ans);

///////////////////// KRIGING /////////////////////////////////////////////////////////////////////////////////////

SEXP kriging(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_CovT, SEXP R_krigType);

/* either kriging all statistics in dual form or with calculating variances
 * and gradient estimation type is the same for all covariance structures  */
SEXP estimateJacobian(SEXP R_Xmat, SEXP R_data, SEXP R_points, SEXP R_covList, SEXP R_krigtype);

SEXP getDualKrigingWeights(SEXP R_Cmat, SEXP R_Fmat, SEXP R_data);

///////////////////// QUASI_LIKELIHOOD  ///////////////////////////////////////////////////////////////////////////

SEXP finalizeQL();

SEXP initQL(SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm);

SEXP quasiDeviance(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat, SEXP R_cm, SEXP R_qdValue);

SEXP mahalanobis(SEXP R_points, SEXP R_qsd, SEXP R_qlopts, SEXP R_X, SEXP R_Vmat,	SEXP R_cm, SEXP R_qdValue);

SEXP qDValue(SEXP R_point);

SEXP mahalValue(SEXP R_point);

//SEXP cvError(SEXP R_point, SEXP R_Xmat, SEXP R_data, SEXP R_covT, SEXP R_krigType, SEXP R_cm);

#endif /* KRIGING_H_ */
