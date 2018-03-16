/**
 * @file        auxil.cc
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions
 *
 *
 */

#include "auxil.h"
#include "error.h"

#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

#define SVD_TOL 1e-12

int check_Lapack_error( int info, const char* name, int line, const char *file);

/** \brief Merge cholesky decomposition into matrix
 *         V=t(L)%*%L
 *
 * @param x  vector of cholesky decomposed entries
 * @param z  result matrix V=t(L)%*%L
 * @param nx dimension
 * @param y  work, nx*nx length
 */

int chol2var(double *x, double *z, int nx, double *y) {
  int info = 0;
  double tmp = 0;
  if(nx>1) {
	  int i,j,k;
      MEMZERO(y,nx*nx);
  	  for(k=j=0; j<nx; j++){
  		  for(i=0;i<j+1; i++,k++){
  			tmp = x[k];
  			y[j*nx+i] = tmp;
  			if (!R_FINITE(tmp))
  			  { info=1; }
  		  }
  	  }
   	  if(info > 0)
   		WRR("`NaN` detected in `chol2var`.");
  	  matmult_trans(y,nx,nx,y,nx,nx,z,info);
  	  if(info > 0)
  		WRR("`NaN` detected in matrix multiplication in `chol2var`.");
  } else {
	  tmp = SQR(*x);
	  if (!R_FINITE(tmp))
		 WRR("`NaN` detected in `chol2var`.");
	  *z = tmp;
  }
  return info;
}


/*! \brief Simple finite difference approximation,
 *    compute gradient or Jacobian
 *
 * @param x vector of coordinates, the point
 * @param nx if x is pointer to a row of a matrix of points, nx are the number of rows of x
 * @param dx length of x
 * @param fval vector of function values at x
 * @param m length of fval
 * @param jac gradient/Jacobian
 * @param fdwork working vector of length m
 * @param func callback function pointer
 * @param data void data pointer
 * @param eps difference
 * @param to_negative whether to multiply by -1
 * @param info integer to signal NA or non finite values
 *
 *
 * @return void
 */
void fdJacobian(double *x, int dx, double *fval, int m, double *jac, double *fdwork,
	  fnCompare_wrap func, void *data, double eps, int to_negative, int &info) {

	int j, k, have_na=0;
	double h, tmp, y;

	for (j=0;j<dx;j++) {
		tmp=x[j];
		h=eps*std::fabs(tmp);
		if(h < DBL_EPSILON)
		 { h=eps; }
		x[j]=tmp+h;
		if (ISNAN(x[j]) || !R_FINITE(x[j]))
		  {have_na = 1; break; }
	    h=x[j]-tmp;										/* finite precision improvement trick */
		func(x,data,fdwork,info);
		if(info > 0)
		  { have_na=1; break;}
		x[j]=tmp;
		if(to_negative)
		 { h = -h; }
		for(k=0;k<m;k++) {
		  y = (fdwork[k] - fval[k])/h;
		  if (ISNAN(y) || !R_FINITE(y))
		    {have_na = 1; break; }
		  //jac[j*m+k] = y;   // original:   m \times n
		  jac[k*dx+j] = y;    // transposed: n rows (parameters) \times m values, e. i. statistics or something else
		}
	}
	info = have_na;
}

void isPositiveDefinite(double *C, int *dim, int *info) {
  int bytes = *dim * *dim;
  double *tmp = CALLOC(bytes, double);
  MEMCPY(tmp, C, bytes);
  F77_CALL(dpofa)(tmp, dim, dim, info);
  FREE(tmp);
}

/*! \brief Distance matrix: (Estimation of Kriging derivatives)
 *              Comment: calculate the distances between
 *                        y_1,...,y_n and all rows x_1,...x_n
 *                        subtract first row of y from all rows of x
 *
 * @param x  matrix, pointer to first element
 * @param nrx rows of x
 * @param xdim columns of x
 * @param y matrix, pointer to some row,
 * @param nry rows of y
 * @param d  distance matrix: size: nrx * ncx
 */
void dist_x0_X( double *x, int nrx, int ncx, double *y, int nry, double *d ) {
	int iy, ix, i, j ,l;

	for(ix=iy=j=l=0; j<ncx; j++, iy+=nry,ix+=nrx ) {
		for(i=0;i<nrx;i++) {
		   d[l++] = std::fabs( *(y+iy) - *(x+ix+i));
		}
	}

}

/** \brief | X1 - X2 | element wise subtraction
 *
 * @param x matrix, pointer to first element matrix, pointer to first element
 * @param nrx rows of x
 * @param xdim columns of x
 * @param y matrix, pointer to some row,
 * @param nry rows of y
 * @param d size: nrx * ncx
 */
void dist_X1_X2( double *x, int nrx, int xdim, double *y, int nry, double *d ) {
	int iy, ix, i, j, l;
	for(ix=iy=j=l=0; j< xdim; j++, iy+=nry,ix+=nrx ) {
		for(i=0;i<nrx;i++) {
			d[l++] = std::fabs( *(y+iy+i) - *(x+ix+i));
		}
	}
}

/**
 * Norm of a vector: (Estimation of Kriging derivatives)
 *  ii: row index of matrix x
 */
double
norm_x_2( double *x, int nrx, int ncx, int ii) {
    int j=0;
    double h=0;
    for(;j<ncx;j++)
       h += SQR(x[ ii + j*nrx ]);
    return std::sqrt(h);
}


/**
 * Distance between two vectors:
 *  x1, x2: matrices
 *  i1,i2: row indices specify which vectors should be used
 */

double norm2( double *x1, int n1, double *x2, int n2, int d, int i1, int i2) {
    double h=0;
    for (int k = 0; k < d; k++){
       h += SQR((x1[i1 + n1 * k] - x2[i2 + n2 * k]));
    }
    return std::sqrt(h);
}

/*
 * x1,x2:    matrices
 * dx1,dx2:  offsets
*/
double norm_2( double *x1, double *x2, int dx1, int dx2, int xdim) {
  int ix, iy, k;
  double h = 0;
  for ( iy = ix = k = 0; k < xdim; k++, iy+=dx2, ix+=dx1 )
    h += SQR(( *(x1+ix) - *(x2+iy)));
  if(!R_FINITE(h)){
	  WRR("`NaN` detected in `norm_2`.")
	  return R_NaN;
  }
  return std::sqrt(h);

}

double denorm(double *x, int n) {
  double value = 0.0;
  for (int i = 0; i < n; ++i )
    value += x[i] * x[i];
  if(!R_FINITE(value)){
  	  WRR("`NaN` detected in `denorm`.")
  	  return R_NaN;
  }
  return std::sqrt(value);
}


double innerProduct(double *x, double *y, int n) {
	double sum=0;
    for(int i=0;i<n;++i)
      sum += x[i]*y[i];
    if(!R_FINITE(sum)) {
      WRR("`NaN` detected in `innerProduct`.")
   	  return R_NaN;
    }
    return sum;
}

/** \brief Set vector and scalar to diagonal of a square matrix
 *
 * @param x    Matrix
 * @param y    Vector
 * @param nx   Dimension
 * @param s    Scalar
 */
void set2diag(double *x, double *y, int *nx, double s) {
  int have_na=0, n=*nx;
  for(int k=0;k<n;++k){
	 x[k*n+k] = y[k] + s;
	 if(!R_FINITE(x[k*n+k]))
	  { have_na=1; break;}
  }
  if(have_na>0)
 	WRR("`NaN` detected in `set2diag`.")
}

int add2diag(double *vmat, int nx, double *s) {
  int have_na=0;
  double *pq = vmat;
  for(int k=0;k<nx; ++k,pq+=nx+1){
	  *pq += s[k];
	  if(!R_FINITE(*pq))
	     { have_na=1; break;}
  }
  if(have_na>0)
 	WRR("`NaN` detected in `add2diag`.")
  return have_na;
}

/**
 * \brief Add Kriging variance of statistics
 *         to diagonal terms of variance matrix
 *
 * @param sig2  kriging variances statistics
 * @param nc    Dimension of vmat
 * @param vmat  Variance matrix
 *
 */
int addVar(double *sig2, int nc, double *vmat, double *work) {
  int have_na=0;
  double *pq = vmat;
  for(int k=0; k<nc; ++k, pq+=nc+1){
    *pq = sig2[k]+work[k];
    if(!R_FINITE(*pq))
      { have_na=1; break;}
  }
  if(have_na>0)
	WRR("`NaN` detected in `addVar`.")
  return have_na;
}

/** \brief Get diagonal of matrix
 *
 * @param x  Matrix
 * @param nx Dimension
 * @param d  Diagonal vector
 */
void splitDiag(double *x, int nx, double *d) {
  for(int k=0;k<nx;k++)
    d[k] = x[k*nx+k];
}

/** \brief Get diagonal value of matrix
 *
 * @param x
 * @param nx
 * @param ii
 * @param ans
 */
void diagValue(double *x, int *nx, int *ii, double *ans) {
  *ans = x[*ii * *nx + *ii];
}

/** \brief Maximum of row sums of matrix
 *
 * @param x
 * @param nrx
 * @param ncx
 * @param sgn   equals either +1 or -1 for diagonal terms
 * @return      maximum of row sums
 */
double matrix_abs_sum(double *x, int *nx, int *sgn) {
  int k=0,l=0, a=*sgn, n=*nx;
  double sum=0, maxSoFar=-HUGE_VAL;

  for(k=0;k<n; ++k) {
    for(l=0,sum=0;l<n;l++) {
      if(l==k)
       sum += a*x[l+n*k];
      else
       sum += std::fabs(x[l+n*k]);
    }
    if(sum>maxSoFar) maxSoFar=sum;
  }
  return maxSoFar;
}


void matrix_col_sum(double *x, int *nx, int *y) {
  int i=0,j=0, n=*nx;
  double sum = 0;
  for(j=0;j<n; ++j) {
      for(y[j]=0,i=0,sum=0;i<n;++i)
         sum += x[i+n*j];
      y[j]=sum;
  }

}

/*! \brief D %*% X
 *
 * @param x X
 * @param pnrx rows X
 * @param pncx col X
 * @param y diagonal matrix stored as a vector
 */
void matmult_diag(double *x, int nrx, int ncx, double *y, int &info) {
	int have_na=0;
	for(int i=0;i<nrx*ncx;i++){
      x[i] *= y[ i%nrx ];
      if(!R_FINITE(x[i]))
        {have_na=1; break; }
    }
	info=have_na;
}

void matmult_diag_sqrt(double *x, int nrx, int ncx, double *y,int &info) {
    int have_na=0;
	for(int i=0;i<nrx*ncx;i++) {
      x[i] *= std::sqrt(y[ i%nrx ]);
      if(!R_FINITE(x[i]))
        {have_na=1; break; }
    }
	info=have_na;
}

void matmult(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z, int &info) {
    int i, j, k,
		have_na = 0;
    long double sum = 0;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
		for (i = 0; i < nrx*ncx; i++)
			if (ISNAN(x[i])) {have_na = 1; break;}
		if (!have_na)
			for (i = 0; i < nry*ncy; i++)
				if (ISNAN(y[i])) {have_na = 1; break;}
		if (!have_na) {
			for (i = 0; i < nrx; i++)
			 for (k = 0; k < ncy; k++) {
				sum = 0.0;
				for (j = 0; j < ncx; j++)
				  sum += x[i + j * nrx] * y[j + k * nry];
				z[i + k * nrx] = sum;
			 }
		} else {
			double one=1.0,
				  zero=0.0;

			F77_CALL(dgemm)("N", "N", &nrx, &ncy, &ncx, &one,
					x, &nrx, y, &nry, &zero, z, &nrx);
		}
	} else /* zero-extent operations should return zeroes */
		for(i = 0; i < nrx*ncy; i++) z[i] = 0;
    info = have_na;
}

/** \brief Matrix multiplication
 *
 * @param x      Matrix, which will to be transposed (not yet transposed at input)
 * @param pnrx   columns of x
 * @param pncx   rows of x
 * @param y      Matrix
 * @param pnry   columns of y
 * @param pncy   columns of y
 * @param z      result matrix
 */
void matmult_trans(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z, int &info)
{
	int i, j, k,
		have_na = 0;
	long double sum = 0;

	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
		for (i = 0; i < nrx*ncx; i++)
		    if (ISNAN(x[i])) {have_na = 1; break;}
		if (!have_na)
		    for (i = 0; i < nry*ncy; i++)
			 if (ISNAN(y[i])) {have_na = 1; break;}
		if (!have_na) {
		    for (i = 0; i < ncx; i++)
			  for (j = 0; j < ncy; j++) {
			    sum = 0.0;
			    for (k = 0; k < nrx; k++)
			    	sum += x[i*nrx+k]*y[j*nrx+k];
			    z[j*ncx+i] = sum;
			  }
		}
    } else /* zero-extent operations should return zeroes */
		for(i = 0; i < nrx*ncy; i++) z[i] = 0;
	info = have_na;
}

int check_Lapack_error(int info, const char* name, int line, const char *file)
{
  if(info == 0)
    return NO_ERROR;
  if(PL > INTERNAL_ERROR)
    Rprintf("error in %s (code=%d) \n %s (line %d).\n", name, info, file, line);
  return LAPACK_ERROR;
}

/*! \brief QR factorisation
 *
 * @param a matrix
 * @param m rows
 * @param n columns
 *
 * @return qr_data object
 */

qr_data qr_dgeqrf(double *a, int m, int n, int &err)
{
    int info=0, lwork = -1;
    double *work=0, tmp=0;

    qr_data qr = CALLOC( 1, struct qr_data_s);
    qr->qr = CALLOC(m*n, double);
    MEMCPY(qr->qr,a,m*n);

    qr->nrow = m;
    qr->ncol = n;
    qr->ntau = ( m < n ) ? m : n;
    qr->tau = CALLOC( qr->ntau, double);

    F77_CALL(dgeqrf)(&(qr->nrow), &(qr->ncol), qr->qr, &(qr->nrow), qr->tau, &tmp, &lwork, &info);
    if( (err=check_Lapack_error(info,"first call to `dgeqrf` failed",__LINE__, __FILE__)) != NO_ERROR )
      goto ErrHandler;

    lwork = (int) tmp;
    work = CALLOC(lwork, double);
    F77_CALL(dgeqrf)(&(qr->nrow), &(qr->ncol), qr->qr, &(qr->nrow), qr->tau, work, &lwork, &info);
    if( (err=check_Lapack_error(info,"second call to `dgeqrf` failed.",__LINE__, __FILE__)) != NO_ERROR )
      goto ErrHandler;

    FREE(work);
    return qr;

ErrHandler:
   FREE(work);
   qrFree(qr);
   err=LAPACK_QR_ERROR;
   LOG_ERROR(LAPACK_QR_ERROR, " `qr_dgeqrf` failed");
   return qr;
}


void qrFree(qr_data qr) {
  if(!qr) return;
  FREE(qr->qr);
  FREE(qr->tau);
  FREE(qr);
}

/*! \brief Compute Projection matrix from QR factorisation
 *
 * @param qr qr_data object
 * @param x (allocated) matrix, on exit the projection matrix
 */

void nullSpaceMat(qr_data qr, double *x, int &err) {
   int i, j, k,
	info = 0, have_na = 0,
      n = qr->ncol,  // cols F
      m = qr->nrow,  // rows F
      n_tau = qr->ntau;
    double *work=0, *Q=0;

    CALLOCX(Q,m*m,double);
    Imat(Q,m);

    int lwork = -1;
    double tmp = 0.0;
    F77_CALL(dormqr)("L", "N", &m, &m, &n_tau, qr->qr, &m, qr->tau, Q, &m, &tmp, &lwork, &info);

    if( (err=check_Lapack_error(info,"first call to `dormqr` failed.",__LINE__, __FILE__)) != NO_ERROR )
       goto ErrHandler;

    lwork = (int) tmp;
    CALLOCX(work,lwork, double);
    F77_CALL(dormqr)("L", "N", &m, &m, &n_tau, qr->qr, &m, qr->tau, Q, &m,  work, &lwork, &info);

    if( (err=check_Lapack_error(info,"second call to `dormqr`",__LINE__, __FILE__)) != NO_ERROR )
       goto ErrHandler;

    //size(x) = m x (m-n)
    //rank(F) has to be maximal!
    for(k=0, j = n; j < m; j++)
      for(i = 0; i < m; i++, k++){
        x[k] = Q[j*m+i];
        if (!R_FINITE(x[i]))
          { have_na = 1; break; }
      }

    err = (have_na > 0 ? NaN_WARNING : NO_ERROR);
	if(have_na > 0)
	   LOG_WARNING(err, " `nullSpaceMat` failed: detected `NaN` in matrix `x`");

    FREE(Q);
    FREE(work);
    return;

ErrHandler:
  FREE(Q);
  FREE(work);
  err=LAPACK_PMAT_ERROR;
  LOG_ERROR(err, " `nullSpaceMat` failed");
}

void solveLU(double *A, int nA, double *B, int nB, int &err) {
  int info = 0;
  int lda = nA, ldb = nA;

  int *ipiv = CALLOC(nA, int);
  F77_NAME(dgesv)(&nA, &nB, A, &lda, ipiv, B, &ldb, &info);
  FREE(ipiv);

  if((err=check_Lapack_error(info,"`dgesv`failed.",__LINE__, __FILE__)) != NO_ERROR){
	  err=LAPACK_SOLVE_ERROR;
	  LOG_ERROR(err, " `solveLU` failed.");
  }

}

void solveQR(double *X, int *nrx, int *ncx, double *y, int *ncy, int &err) {
  int info = 0, lwork=-1;
  int  nrowx = *nrx, ncolx = *ncx, ncoly = *ncy;
  double *work=0, *Xtmp=0, wtmp=0;

  CALLOCX(Xtmp,nrowx*ncolx, double);
  MEMCPY(Xtmp,X,nrowx*ncolx);

  F77_CALL(dgels)("N", &nrowx, &ncolx, &ncoly, Xtmp,
                  &nrowx, y, &nrowx, &wtmp, &lwork, &info);

  if( (err=check_Lapack_error(info," first call to `dgels failed!",__LINE__, __FILE__)) != NO_ERROR)
     goto ErrHandler;

  lwork = (int) wtmp;
  CALLOCX(work,lwork, double);
  F77_CALL(dgels)("N", &nrowx, &ncolx, &ncoly, Xtmp,
                  &nrowx, y, &nrowx, work, &lwork, &info);

  if( (err=check_Lapack_error(info,"second call to `dgels` failed!",__LINE__, __FILE__)) != NO_ERROR)
    goto ErrHandler;

  FREE(work);
  FREE(Xtmp);
  return;

ErrHandler:
   FREE(Xtmp);
   FREE(work);
   err=LAPACK_QR_ERROR;
   LOG_ERROR(err, " `solveQR` failed");
}

/**
 *  \brief Cholesky solve
 *   Comment: A is overwritten here!
 */
void solveCH(double *A, int n, double *y, int ncoly, int &err) {
    int info=0;
    /* factorization */
	F77_NAME(dpotrf)("U",&n,A,&n,&info);
	if((err=check_Lapack_error(info," Cholesky factorization `dpotrf` failed.", __LINE__,__FILE__)) != NO_ERROR){
		err=LAPACK_FACTORIZE_ERROR;
		LOG_ERROR(err, " `solveCH` failed: matrix is probably not positive (semi)definite.");
		return;
	}

	F77_NAME(dpotrs)("U", &n, &ncoly, A, &n, y, &n, &info);
	if((err=check_Lapack_error(info," `dpotrs` failed!",__LINE__, __FILE__)) != NO_ERROR){
		err=LAPACK_SOLVE_ERROR;
		LOG_ERROR(err, " `solveCH` failed: solving for solution matrix failed.");
	}
}


void solveDSPTRS(double *A, int n, double *B, int nrhs, int &err) {
  int info = 0,
	  nAP = n*(n+1)/2;
  double *AP = CALLOC(nAP, double);

  /* Maybe we make use of packed storage pattern later? */
  triangMat_U(A,n,AP, info);
  if(info > 0){
    FREE(AP);
    err = NaN_WARNING;
  	LOG_WARNING(err," `solveDSPTRS` failed: `NaNs` detected in `triangMat_U`.")
    return;
  }

  /* Matrix factorization Bunch-Kaufmann */
  int *ipiv = CALLOC(n,int);
  F77_NAME(dsptrf)("U",&n,AP,ipiv,&info);
  if((err=check_Lapack_error(info," `dsptrf` failed!",__LINE__, __FILE__)) != NO_ERROR){
	  err=LAPACK_FACTORIZE_ERROR;
	  LOG_ERROR(err, " `solveDSPTRS` failed");
	  goto ErrHandler;
  }

  /* Solving the linear equation */
  F77_NAME(dsptrs)("U",&n,&nrhs,AP,ipiv,B,&n,&info);
  if((err=check_Lapack_error(info," `dsptrs` failed!",__LINE__, __FILE__)) != NO_ERROR){
	  err=LAPACK_SOLVE_ERROR;
	  LOG_ERROR(err, " `solveDSPTRS` failed");
	  goto ErrHandler;
  }

  FREE(AP);
  FREE(ipiv);
  return;

ErrHandler:
  FREE(AP);
  FREE(ipiv);

}


void gsiSolve(double *A, int n, double *B, int nrhs, double *Awork,
				int &err, inversion_type type) {

	int i=0, have_na=0,
		info=0, n2=SQR(n);

	if(type == Chol){
		// copy input matrix
		for(i=0; i < n2; i++) {
		  if(!R_FINITE(A[i])) { have_na=1; break; }
		  Awork[i] = A[i];
		}
		err = (have_na > 0 ? NaN_WARNING : NO_ERROR);
		if(have_na > 0){
		  LOG_WARNING(err, " `gsiSolve` failed  by `dpotrs`: detected `NaN` in matrix `A`.");
		  return;
		}

		// factorization
		F77_NAME(dpotrf)("U",&n,Awork,&n,&info);
		if((err=check_Lapack_error(info," `dpotrf` failed.", __LINE__,__FILE__)) == NO_ERROR) {
			F77_NAME(dpotrs)("U", &n, &nrhs, Awork, &n, B, &n, &info);
			err=check_Lapack_error(info," `dpotrs` failed.", __LINE__,__FILE__);
		}
		// try `SVD` if `Chol` has failed
		if(err != NO_ERROR) {
			type = SVD;
			err = LAPACK_SOLVE_ERROR;
			if(PL > INTERNAL_ERROR)
			  LOG_WARNING(err, " `gsiSolve` failed by Cholesky factorization. Matrix is probably not positive (semi)definite.");
		}
	}

	if(type == SVD)	{
		double wkopt=0, *work=&wkopt,
			  *s=0, *u=0, *vt=0, *vwork=0;
		int *iwork=0, size8=n*8;

		// copy input matrix
		have_na=0;
		for(i=0; i < n2; i++) {
		  if(!R_FINITE(A[i]))
		    { have_na=1; break; }
		  Awork[i] = A[i];
		}
		err = (have_na > 0 ? NaN_WARNING : NO_ERROR);
		if(have_na > 0){
		   LOG_WARNING(err, " `gsiSolve` failed: detected `NaN` in matrix `A`");
		   return;
		}

		CALLOCX(s,n,double);
		CALLOCX(u,n2,double);
		CALLOCX(vt,n2,double);
		CALLOCX(iwork, size8, int);

		int lwork = -1;
		F77_CALL(dgesdd)("A", &n, &n, Awork, &n, s, u, &n, vt, &n, work, &lwork, iwork, &info);
		if((err=check_Lapack_error(info," first call to `dgesdd` failed.",__LINE__, __FILE__)) == NO_ERROR ){
			lwork = (int) wkopt;
			CALLOCX(work,lwork,double);

			F77_CALL(dgesdd)("A", &n, &n, Awork, &n, s, u, &n, vt, &n, work, &lwork, iwork, &info);
			if((err=check_Lapack_error(info," second call to `dgesdd` failed.",__LINE__, __FILE__)) == NO_ERROR ){
				CALLOCX(vwork,n*nrhs,double);
				matmult(vt,n,n,B,n,nrhs,vwork,info);
				if(info > 0){
					err=NaN_WARNING;
					WRR("`NaN` detected in matrix multiplication.")
				}

				/* invert diagonal terms in matrix vwork */
				for(int i=0; i < n; i++){
				  if(s[i] < SVD_TOL) {
					for(int j=0; j < nrhs; j++)
					 vwork[j*n+i] = 0.0;
				  } else {
					for(int j=0; j < nrhs; j++) {
					  if(!R_FINITE(s[i])) { have_na=1; break; }
					  vwork[j*n+i] /= s[i];
					}
				  }
				}

				err = (have_na > 0 ? NaN_WARNING : NO_ERROR);
				if(have_na > 0)
					LOG_WARNING(err," `gsiSolve` failed by `dgesdd`: detected `NaN` in matrix `vwork`.");

				/** and multiply from right by D*U */
				matmult(u,n,n,vwork,n,nrhs,B,info);
				if(info > 0){
				  err=NaN_WARNING;
			      WRR("`NaN` detected in matrix multiplication.")
				}
				FREE(vwork);
			}
			FREE(work);
		}
		FREE(s);
		FREE(u);
		FREE(vt);
		FREE(iwork);

		if(err != NO_ERROR){
			type = Bunch;
			err = LAPACK_SOLVE_ERROR;
			if(PL > INTERNAL_ERROR)
			  LOG_ERROR(err, " `gsiSolve` failed by SVD.");
		}
	}

	if(type == Bunch) {
		err = NO_ERROR;
		solveDSPTRS(A,n,B,nrhs,err);
		if(err != NO_ERROR){
			err = LAPACK_SOLVE_ERROR;
			if(PL > 0)
			 LOG_ERROR(err, " `gsiSolve` finally failed by Bunch.");
		}
	}

}

///** \brief Factorize matrix A as A=LL^T
// *
// * @param A     matrix
// * @param nA    rows of A
// * @param info  error info from subsequent call
// */
//void factorize_chol_L(double *A, int *nA, int *err) {
//  *err=NO_ERROR;
//  F77_NAME(dpotrf)("L",nA,A,nA,err);
//}

void factorize_chol_L(double *A, int *nA, int &err) {
   int info = 0, n = *nA;
   const char *uplo = "L";

   /* factorize */
   F77_NAME(dpotrf)(uplo,&n,A,&n,&info);
   if((err=check_Lapack_error(info,"call to 'dpotrf' failed!",__LINE__, __FILE__)) != NO_ERROR ){
	   err=LAPACK_FACTORIZE_ERROR;
	   LOG_ERROR(err, "'factorize_chol_L' failed");
   }
}

/** \brief Solve AX=B with A=LL^T
 *         Only the lower triangular of A is used
 *
 * @param A     matrix A, stored as array
 * @param nA    rows of A
 * @param B     matrix B, stored as array
 * @param nB    columns of B, number of right hand sides (nrhs)
 * @param err   error indicator
 */
void solve_chol_factorized(double *A, int *nA, double *B, int *nB, int &err) {
  int info = 0,  n = *nA, nrhs = *nB;  // B has dimension:  n x nrhs, where nrhs= #cols
  const char *uplo = "L";

  /* solve with factorized matrix A=LL^T */
  F77_NAME(dpotrs)(uplo,&n,&nrhs,A,&n,B,&n,&info);
  if((err=check_Lapack_error(info,"call to 'dpotrs' failed! ",__LINE__, __FILE__)) != NO_ERROR ){
	  err=LAPACK_SOLVE_ERROR;
	  LOG_ERROR(err, " 'solve_chol_factorized' failed!");
  }

}


void solve_chol_triangular(double *A, int *nA, double *B, int *nB, int &err) {
  int info = 0,  n = *nA, nrhs = *nB;  // B has dimension:  n x nrhs, where nrhs= #cols
  const char *diag = "N";
  const char *tran = "N";
  const char *uplo = "L";

  F77_NAME(dtrtrs)(uplo,tran,diag,&n,&nrhs,A,&n,B,&n,&info);
  if((err=check_Lapack_error(info,"call to 'dtrtrs' failed! ",__LINE__, __FILE__)) != NO_ERROR ){
	  err=LAPACK_SOLVE_ERROR;
	  LOG_ERROR(err, " 'solve_chol_factorized' failed!");
  }
}

void invMatrix(double *A, int n, double *ans, int &err, inversion_type type) {
	const char *uplo="U";
	int info=0, have_na=0, n2=SQR(n);

	if(type == Chol){
		// copy input matrix
		MEMCPY(ans,A,n2);

		// factorization
		F77_NAME(dpotrf)(uplo,&n,ans,&n,&info);
		if((err=check_Lapack_error(info," `dpotrf` failed.", __LINE__,__FILE__)) == NO_ERROR) {
			/* invert */
			F77_NAME(dpotri)(uplo, &n, ans, &n, &info);
			if((err=check_Lapack_error(info," `dpotri` failed.", __LINE__,__FILE__)) == NO_ERROR){
				int i2, i3, j, i;
				for (i2 = i = 0; i<n; i++, i2 += n+1) {
				  for (i3 = i2+1, j = i2+n; j < n2; j += n){
					  if (!R_FINITE(ans[i]))
						{ have_na = 1; break; }
					  ans[i3++] = ans[j];
				  }
				}
				err = (have_na > 0 ? NaN_WARNING : NO_ERROR);
				if(have_na > 0)
				  LOG_WARNING(err," `invMatrix` failed.");
			}
		}
		/// try `SVD` if `Chol` has failed
		if(err != NO_ERROR) {
			type = SVD;
			err = LAPACK_INVERSION_ERROR;
			if(PL > INTERNAL_ERROR)
			  LOG_ERROR(err, " `invMatrix` failed by Cholesky factorization. Matrix is probably not positive (semi)definite.");
		}
	}

	if(type == SVD)	{
		double wkopt=0, *work=&wkopt,
			  *s=0, *u=0, *vt=0;
		int *iwork=0, size8=n*8;
		have_na = 0;
		// copy input matrix
		MEMCPY(ans,A,n2);

		CALLOCX(s,n,double);
		CALLOCX(u,n2,double);
		CALLOCX(vt,n2,double);
		CALLOCX(iwork, size8, int);

		int lwork = -1;
		F77_CALL(dgesdd)("A", &n, &n, ans, &n, s, u, &n, vt, &n, work, &lwork, iwork, &info);
		if((err=check_Lapack_error(info," first call to `dgesdd` failed.",__LINE__, __FILE__)) == NO_ERROR ){
			lwork = (int) wkopt;
			CALLOCX(work,lwork,double);

			F77_CALL(dgesdd)("A", &n, &n, ans, &n, s, u, &n, vt, &n, work, &lwork, iwork, &info);
			if((err=check_Lapack_error(info," second call to `dgesdd` failed.",__LINE__, __FILE__)) == NO_ERROR ){
				/* invert diagonal terms in matrix Vt */
				for(int i=0; i < n; i++){
				  if(s[i] < SVD_TOL) {
					for(int j=0; j < n; j++)
					 vt[j*n+i] = 0.0;
				  } else {
					for(int j=0; j < n; j++) {
					  if(!R_FINITE(s[i])) { have_na=1; break; }
					  vt[j*n+i] /= s[i];
					}
				  }
				}
				err = (have_na > 0 ? NaN_WARNING : NO_ERROR);
				if(have_na > 0)
				  LOG_WARNING(err," in `invMatrix`.");

				/** and multiply from right by D*U */
				matmult(u,n,n,vt,n,n,ans,info);
				if(info > 0){
					err=NaN_WARNING;
					WRR("`NaN` detected in matrix multiplication.")
				}
			}
			FREE(work);
		}
		FREE(s);
		FREE(u);
		FREE(vt);
		FREE(iwork);

		if(err != NO_ERROR){
			type = Bunch;
			err = LAPACK_INVERSION_ERROR;
			if(PL > INTERNAL_ERROR)
				LOG_ERROR(err, " `invMatrix` failed by SVD.");
		}
	}

	if(type == Bunch)
	{
		int nAP = n*(n+1)/2;
		double *AP = NULL;
		CALLOCX(AP,nAP,double);

		have_na = 0;
		triangMat_U(A,n,AP,info);
		if(info > 0){
		  FREE(AP);
		  err = NaN_WARNING;
		  LOG_WARNING(err,"`invMatrix` failed: `NaNs` detected in `triangMat_U`.")
		}

		int *ipiv = NULL;
		CALLOCX(ipiv,n,int);

		/* Matrix factorization */
		F77_NAME(dsptrf)(uplo,&n,AP,ipiv,&info);
		if((err=check_Lapack_error(info," matrix factorization `dsptrf` failed.",__LINE__, __FILE__)) == NO_ERROR ){
			double *work = NULL;
			CALLOCX(work,n,double);

			/* Inverse, possibly indefinite */
			F77_NAME(dsptri)(uplo,&n,AP,ipiv,work,&info);
			if((err=check_Lapack_error(info," Inversion by `dsptri` failed!", __LINE__,__FILE__)) == NO_ERROR){
				triangMat_U_back(ans,n,AP,info); // result in ans
				if(info > 0){
					err =  NaN_WARNING;
					LOG_WARNING(err," `invMatrix` failed: `NaNs` detected in `triangMat_U_back`.");
				}
			}
			FREE(work);
		}

		if(err != NO_ERROR) {
			err = LAPACK_INVERSION_ERROR;
			LOG_ERROR(err, " `invMatrix` failed by Bunch-Kaufman factorization.");
		}

		FREE(AP);
		FREE(ipiv);
	}

}

/**
 * @brief: Transpose matrix: Columns of x to rows of y
 *
 * @param y
 * @param dy
 * @param x
 * @param dx
 * @param nr
 * @param nc
 */
void mat_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol, int &info){
  for (int i=0; i<nrow; i++) {
    for (int j=0;j<ncol; j++) {
       y[j] = x[i + j * ldx]; // column of x to row of y
       if (ISNAN(y[i]))
        { info = 1; break; }
    }
    y += ldy;
  }
}

void triangMat_U(double *A, int nA, double *AP, int &info) {
    int n1 = nA, i, j, k;

    for(info=k=j=0; j<n1; j++)
    	for(i=0;i<j+1; i++,k++){
    		if (ISNAN(A[i]))
    		  { info = 1; break; }
    		AP[k] = A[j*n1+i];
    	}
}

void triangMat_U_back(double *A, int nA, double *AP, int &info) {
    int n1 = nA, i, j, k;

    for(info=k=j=0; j<n1; j++)
    	for(i=0;i<j+1; i++,k++){
    		if (ISNAN(AP[i]))
    		 { info = 1; break; }
    		A[i*n1+j] = A[j*n1+i] = AP[k];
    	}
}

void Imat(double *x, int n) {
  // matrix must be initialized by zero!!
  int i;
  for(i=0;i<n*n;i++)
    x[i]=0.0;
  for(i=0;i<n;i++)
    x[i+n*i]=1.0;
}

void trendfunc(double *x, int nx, int dx, double *f, int model){
	if(!model)
	 Rf_error(_("No trend order specified!\n"));

	switch (model) {
	case 0:
		basis0(f);
		break;
	case 1:
		basis1(x,nx,dx,f);
		break;
	case 2:
		basis2(x,nx,dx,f);
		break;
	default:
		ERR("Invalid trend order number. Choose either k=1,2.");
	    break;
	};

}

void basis0(double *f) {
	f[0] = 1;
}

void basis1(double *x, int inx, int dx, double *f)
{
	f[0] = 1;
	int i, ix;
	for(ix = i = 0; i < dx; i++, ix+=inx) f[i+1] = *(x+ix);
}

void basis2( double *x, int inx, int dx, double *f)
{
	int k, i, ix, iy = 0, l = dx+1;

	f[0] = 1;
	for(ix = i = 0; i < dx; i++, ix+=inx) f[i+1] = *(x+ix);

	for (ix = k = 0; k < dx; k++, ix+=inx)
	   for ( iy = ix, i = k; i < dx; i++ , iy+=inx)
		 f[l++] = *(x+ix) * *(x+iy);
}

void Fmatrix(double *x,double *F, int n, int d, int trend) {
	/** Size F matrix: stored columns first
	 *
	 * constant trend:  size(x) x 1
	 * linear trend  :  size(x) x (dim + 1)
	 * quadratic trend: size(x) x (dim+1)*(dim+2)/2
	 *
	 */

	double *pF = F;

	if(trend >= 0 && trend < 3)
	  //constant terms
	    if(n < 2)
	      Rf_error(_("Trend matrix is singular!\n"));
	for(int i = 0; i < n; i++, pF++) *pF = 1;
	if(trend > 0) {
		//linear terms
	   if(n < d+2)
	     Rf_error(_("Trend matrix is singular!\n"));

	   for(int j = 0; j < d; j++) //columns of x
		   for(int k = 0; k < n; k++,pF++) // rows of F
			    *pF = x[j*n+k];
	}
	if(trend > 1) {
	 //quadratic terms
		if( n < (d+1)*(d+2)/2+1 )
	      Rf_error(_("Trend matrix is singular!\n"));

		for (int k = 0; k < d; k++)
		  for (int l = k; l < d; l++)
		    for ( int i = 0; i < n; i++, pF++)
			  *pF = x[k*n+i]*x[l*n+i];

	} else if(trend > 2){
	    ERR("Invalid trend order number. Choose either k=1,2.");
	}

}
