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

int check_Lapack_error( int info, const char* name, int line, const char *file);

/** \brief Merge cholesky decomposition into matrix
 *         V=t(L)%*%L
 *
 * @param x  vector of cholesky decomposed entries
 * @param z  result matrix V=t(L)%*%L
 * @param nx dimension
 * @param y  work, nx*nx length
 */

void chol2var(double *x, double *z, int nx, double *y) {
  int info = 0;
  double tmp = 0;
  if(nx>1) {
	  int i,j,k;
      MEMZERO(y,nx*nx);
  	  for(k=j=0; j<nx; j++)
  		  for(i=0;i<j+1; i++,k++){
  			  tmp = x[k];
  			  if (!R_finite(tmp) || ISNA(tmp) || ISNAN(tmp))
  			    WRR("`NaN` detected.");
  			  y[j*nx+i] = tmp;
  		  }
  	  matmult_trans(y,nx,nx,y,nx,nx,z,&info);
  	  if(info)
  		WRR("`NaN` detected in matrix multiplication.");
  } else {
	  tmp = SQR(*x);
	  if (!R_finite(tmp) || ISNA(tmp) || ISNAN(tmp))
		  WRR("`NaN` detected.");
	  *z = tmp;
  }
}


/*! \brief Finite difference approximation,
 *            Comment: Calculate gradient or jacobian
 *
 * @param x vector of coordinates, the point
 * @param nx if x is pointer to a row of a matrix of points, nx are the number of rows of x
 * @param dx length of x
 * @param fval vector of function values at x
 * @param m length of fval
 * @param jac gradient/jacobian
 * @param func callback function pointer
 * @param data void data pointer
 * @param eps difference
 *
 *
 * @return void
 */
void fdJac(double *x, int dx, double *fval, int m, double *jac,
            fnCompare func, void *data, double eps, int to_negative, int *info) {

	int j, k, have_na = 0;
	double h, temp, *fval_tmp, *x_tmp, y;

	CALLOCX(fval_tmp,m,double);
	CALLOCX(x_tmp,dx,double);
	MEMCPY(x_tmp,x,dx);

	for (j=0;j<dx;j++) {
		temp=x_tmp[j];
		h=eps*fabs(temp);
		if(h==0) h=eps;
		 x_tmp[j]=temp+h;
		 if (ISNAN(x_tmp[j]) || !R_FINITE(x_tmp[j]))
		   {have_na = 1; break; }
		 h=x_tmp[j]-temp;
		 func(x_tmp, data, fval_tmp);
		 x_tmp[j]=temp;
		 if(to_negative) h = -h;
		 for(k=0;k<m;k++) {
		  y = (fval_tmp[k] - fval[k])/h;
		  if (ISNAN(y) || !R_FINITE(y))
		    {have_na = 1; break; }
		  jac[k*dx+j] = y;
		 }
	}
	FREE(x_tmp);
	FREE(fval_tmp);
	*info = have_na;
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
		   d[l++] = fabs( *(y+iy) - *(x+ix+i));
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
			d[l++] = fabs( *(y+iy+i) - *(x+ix+i));
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
  return std::sqrt(h);

}

double denorm (double *x, int n) {
  double value = 0.0;
  for (int i = 0; i < n; i++ )
    value += x[i] * x[i];
  return std::sqrt(value);
}


double innerProduct(double *x, double *y, int n) {
    int i=0;
    double sum=0;
    for(;i<n;++i)
      sum += x[i]*y[i];
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
  int k=0, n=*nx;
  for(;k<n;++k)
    x[k*n+k] = y[k] + s;
}

void add2diag(double *vmat, int nx, double *s) {
  int k=0;
  double *pq = vmat;
  for(;k<nx; ++k,pq+=nx+1)
    *pq += s[k];
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
void
addVar(double *sig2, int nc, double *vmat, double *work) {
  int k=0;
  double *pq = vmat;
  for(;k<nc; ++k, pq+=nc+1)
    *pq = sig2[k]+work[k];
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
       sum += fabs(x[l+n*k]);
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
void matmult_diag(double *x, int nrx, int ncx, double *y) {
    for(int i=0;i<nrx*ncx;i++)
      x[i] *= y[ i%nrx ];
}

void matmult_diag_sqrt(double *x, int nrx, int ncx, double *y) {
    for(int i=0;i<nrx*ncx;i++) {
      x[i] *= std::sqrt(y[ i%nrx ]);
    }
}

void matmult(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z, int *info) {
    const char *transa = "N", *transb = "N";
    int i, j, k,
		have_na = 0;
    double one = 1.0, zero = 0.0;
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
		} else
			F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
					x, &nrx, y, &nry, &zero, z, &nrx);
	} else /* zero-extent operations should return zeroes */
		for(i = 0; i < nrx*ncy; i++) z[i] = 0;
    *info = have_na;
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
void matmult_trans(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z, int *info)
{
	int i, j, k,
		have_na = 0;

	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
		for (i = 0; i < nrx*ncx; i++)
		    if (ISNAN(x[i])) {have_na = 1; break;}
		if (!have_na)
		    for (i = 0; i < nry*ncy; i++)
			if (ISNAN(y[i])) {have_na = 1; break;}
		if (!have_na) {
		    for (i = 0; i < ncx; i++)
			  for (j = 0; j < ncy; j++) {
				z[j*ncx+i] = 0;
			    for (k = 0; k < nrx; k++)
			    	z[j*ncx+i] += x[i*nrx+k]*y[j*nrx+k];
			  }
		}
    } else /* zero-extent operations should return zeroes */
		for(i = 0; i < nrx*ncy; i++) z[i] = 0;
	*info = have_na;
}

int check_Lapack_error(int info, const char* name, int line, const char *file)
{
  if (info == 0)
    return NO_ERROR;

  char MSG[100]="";
  if (info < 0)
    std::sprintf(MSG, "Argument %d to Lapack function %s is illegal. \n ", -info, name);
  else if(info > 0)
    std::sprintf(MSG, "Lapack function %s returned error code %d \n ", name, info);

  PRINT_MSG(MSG);

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

qr_data qr_dgeqrf(double *a, int m, int n, int *err)
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

    if( (*err=check_Lapack_error(info,"First call to dgeqrf",__LINE__, __FILE__)) != NO_ERROR )
      goto ErrHandler;

    lwork = (int) tmp;
    work = CALLOC(lwork, double);
    F77_CALL(dgeqrf)(&(qr->nrow), &(qr->ncol), qr->qr, &(qr->nrow), qr->tau, work, &lwork, &info);

    if( (*err=check_Lapack_error(info,"Second call to dgeqrf",__LINE__, __FILE__)) != NO_ERROR )
      goto ErrHandler;

    FREE(work);
    return qr;

ErrHandler:
   qrFree(qr);
   FREE(work);

   std::strcpy(ERROR_LOC, __FILE__);
   *err=LAPACK_QR_ERROR;
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

void nullSpaceMat(qr_data qr, double *x, int *err) {
    int i,j,k,info;
    int n = qr->ncol;  // cols F
    int m = qr->nrow;  // rows F
    int n_tau = qr->ntau;
    double *work, *Q;

    Q = CALLOC(m*m,double);
    Imat(Q,m);

    int lwork = -1;
    double tmp = 0.0;
    F77_CALL(dormqr)("L", "N", &m, &m, &n_tau, qr->qr, &m, qr->tau, Q, &m, &tmp, &lwork, &info);

    if( (*err=check_Lapack_error(info,"First call to dormqr",__LINE__, __FILE__)) != NO_ERROR )
       goto ErrHandler;

    lwork = (int) tmp;
    work = CALLOC(lwork, double);
    F77_CALL(dormqr)("L", "N", &m, &m, &n_tau, qr->qr, &m, qr->tau, Q, &m,  work, &lwork, &info);
    if(work) FREE(work);

    if( (*err=check_Lapack_error(info,"Second call to dormqr",__LINE__, __FILE__)) != NO_ERROR )
       goto ErrHandler;

    //size(x) = m x (m-n)
    //rank(F) has to be maximal!
    for(k=0, j = n; j < m;j++)
      for(i = 0; i < m; i++, k++)
        x[k] = Q[j*m+i];

    FREE(Q);
    return;

ErrHandler:
  FREE(Q);
  FREE(work);
  std::strcpy(ERROR_LOC, __FILE__);
  *err=LAPACK_PMAT_ERROR;
}


void solveLU(double *A, int nA, double *B, int nB, int *err) {
  int info = 0;
  int lda = nA, ldb = nA;

  int *ipiv = CALLOC(nA, int);
  F77_NAME(dgesv)(&nA, &nB, A, &lda, ipiv, B, &ldb, &info);
  FREE(ipiv);

  if((*err=check_Lapack_error(info,"LU solve failed in dgesv!",__LINE__, __FILE__)) != NO_ERROR) {
      std::strcpy(ERROR_LOC, __FILE__);
      *err = LAPACK_SOLVE_ERROR;
  }
}

void solveQR(double *X, int *nrx, int *ncx, double *y, int *ncy, int *err) {
  int info = 0, lwork=-1;
  int  nrowx = *nrx, ncolx = *ncx, ncoly = *ncy;
  double *work=0,*xtmp=0, tmp=0, *ans=0;

  const char transx = 'N';
  xtmp = CALLOC(nrowx*ncolx, double);
  MEMCPY(xtmp,X,nrowx*ncolx);

  ans = y;
  F77_CALL(dgels)(&transx, &nrowx, &ncolx, &ncoly, xtmp,
                  &nrowx, ans, &nrowx, &tmp, &lwork, &info);

  if( (*err=check_Lapack_error(info,"First call to dgels in solveQR failed!",__LINE__, __FILE__)) != NO_ERROR)
     goto ErrHandler;

  lwork = (int) tmp;
  work = CALLOC(lwork, double);
  F77_CALL(dgels)(&transx, &nrowx, &ncolx, &ncoly, xtmp,
                  &nrowx, ans, &nrowx, work, &lwork, &info);

  if( (*err=check_Lapack_error(info,"Second call to dgels in solveQR failed!",__LINE__, __FILE__)) != NO_ERROR)
    goto ErrHandler;

  FREE(work);
  FREE(xtmp);
  return;

ErrHandler:
   FREE(xtmp);
   FREE(work);
   std::strcpy(ERROR_LOC, __FILE__);
   *err = LAPACK_SOLVE_ERROR;
}

/*! cholesky solve */
void solveCH(double *X, int nrowx, int ncolx, double *y, int ncoly, double *ans /* ncx x ncy */, int *err) {
    int  info = 0;
    const double zero = .0, one = 1.0;

    F77_CALL(dgemm)("T", "N", &nrowx, &ncoly, &nrowx, &one, X, &nrowx, y, &nrowx,&zero, ans, &ncolx);

    double *xtmp = CALLOC(ncolx*ncolx, double);
    F77_CALL(dsyrk)("U", "T", &ncolx, &nrowx, &one, X, &nrowx, &zero,  xtmp, &ncolx);
    F77_CALL(dposv)("U", &ncolx, &ncoly, xtmp, &ncolx, ans, &ncolx, &info);
    FREE(xtmp);

    if( (*err=check_Lapack_error(info," `dposv` failed!",__LINE__, __FILE__)) != NO_ERROR) {
       LOG_ERROR(LAPACK_SOLVE_ERROR, info, " `solveCH` failed");
    }
}

/*! packed storage solve */
void solve_DSPTRS(double *A, int n, double *B, int nrhs, int *err ) {
  *err=NO_ERROR;
  const char *uplo = "U";

  int nAP = n*(n+1)/2;
  int *ipiv = CALLOC(n,int);
  double *AP = CALLOC(nAP, double);

  /* Maybe we make use of packed storage pattern later? */
  triangMat_U(A,n,AP,nAP);

  int info = 0;
  /* Matrix factorization */
  F77_NAME(dsptrf)(uplo,&n,AP,ipiv,&info);
  if(check_Lapack_error(info,"First call to dsptrf in solve_DSPTRS failed!",__LINE__, __FILE__) != NO_ERROR)
    goto ErrHandler;

  /* Solving the linear equation */
  F77_NAME(dsptrs)(uplo,&n,&nrhs,AP,ipiv,B,&n,&info);
  if(check_Lapack_error(info,"Second call to dsptrf in solve_DSPTRS failed!",__LINE__, __FILE__) != NO_ERROR)
      goto ErrHandler;

  FREE(ipiv);
  FREE(AP);
  return;

ErrHandler:
    FREE(AP);
    FREE(ipiv);

  LOG_ERROR(LAPACK_SOLVE_ERROR, info,"Routine solve_DSPTRS failed");
  *err=LAPACK_SOLVE_ERROR;
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

void factorize_chol_L(double *A, int *nA, int *err) {
   int info = 0, n = *nA;
   const char *uplo = "L";
   *err=NO_ERROR;

   /* factorize */
   F77_NAME(dpotrf)(uplo,&n,A,&n,&info);

   if(check_Lapack_error(info,"call to 'dpotrf' failed!",__LINE__, __FILE__) != NO_ERROR )
      goto ErrHandler;

   return;

ErrHandler:
 LOG_ERROR(LAPACK_FACTORIZE_ERROR, info,"'factorize_chol_L' failed");
 *err=LAPACK_FACTORIZE_ERROR;
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
void solve_chol_factorized(double *A, int *nA, double *B, int *nB, int *err) {
  int info = 0,  n = *nA, nrhs = *nB;  // B has dimension:  n x nrhs, where nrhs= #cols
  const char *uplo = "L";
  *err=NO_ERROR;

  /* solve with factorized matrix A=LL^T */
  F77_NAME(dpotrs)(uplo,&n,&nrhs,A,&n,B,&n,&info);
  if(check_Lapack_error(info,"call to 'dpotrs' failed! ",__LINE__, __FILE__) != NO_ERROR )
    goto ErrHandler;

  return;

ErrHandler:
  LOG_ERROR(LAPACK_SOLVE_ERROR, info," 'solve_chol_factorized' failed!");
  *err=LAPACK_SOLVE_ERROR;
}


void solve_chol_triangular(double *A, int *nA, double *B, int *nB, int *err) {
  int info = 0,  n = *nA, nrhs = *nB;  // B has dimension:  n x nrhs, where nrhs= #cols
  const char *diag = "N";
  const char *tran = "N";
  const char *uplo = "L";
  *err=NO_ERROR;

  F77_NAME(dtrtrs)(uplo,tran,diag,&n,&nrhs,A,&n,B,&n,&info);

  if(check_Lapack_error(info,"call to 'dtrtrs' failed! ",__LINE__, __FILE__) != NO_ERROR )
     goto ErrHandler;

   return;

 ErrHandler:
   LOG_ERROR(LAPACK_SOLVE_ERROR, info," 'solve_chol_factorized' failed!");
   *err=LAPACK_SOLVE_ERROR;
}


void invMatrix(double *A,int nA, int *err) {
	int info = 0;
	const char *uplo = "U";
	*err=NO_ERROR;

	int nAP = nA*(nA+1)/2;
	double *AP,*work;
	CALLOCX(AP,nAP,double);
	CALLOCX(work,nA, double);

	triangMat_U(A,nA,AP,nAP);

	int *ipiv;
	CALLOCX(ipiv,nA, int);

	/* Matrix factorization */
	F77_NAME(dsptrf)(uplo,&nA,AP,ipiv,&info);
	if(check_Lapack_error(info,"Matrix factorization dsptrf failed!",__LINE__, __FILE__) != NO_ERROR )
	    goto ErrHandler;

	/* Inverse, possibly indefinite */
	F77_NAME(dsptri)(uplo,&nA,AP,ipiv,work,&info);
	if(check_Lapack_error(info,"Matrix inversion dsptri failed!", __LINE__,__FILE__) != NO_ERROR)
	    goto ErrHandler;

	triangMat_U_back(A,nA,AP,nAP);

	FREE(ipiv);
	FREE(AP);
	FREE(work);

	return;

ErrHandler:
  FREE(ipiv);
  FREE(AP);
  FREE(work);

  LOG_ERROR(LAPACK_INVERSION_ERROR, info,"Routine invMatrix failed");
  *err=LAPACK_INVERSION_ERROR;

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
void mat_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol)
{
  for (int i=0; i<nrow; i++) {
    for (int j=0;j<ncol; j++)
        y[j] = x[i + j * ldx]; // column of x to row of y
    y += ldy;
  }
}

void triangMat_U(double *A, int nA, double *AP, int nAP) {
    int n1 = nA, n2 = nAP, i, j, k;

    if( (n1*(n1+1)/2) != n2 )
    	error(_("Invalid array dimension!\n"));

    for(k=j=0; j<n1; j++)
    	for(i=0;i<j+1; i++,k++)
    	   AP[k] = A[j*n1+i];
}

void triangMat_U_back(double *A, int nA, double *AP, int nAP) {
    int n1 = nA, n2 = nAP, i, j, k;
    if( (n1*(n1+1)/2) != n2 )
    	error(_("Invalid array dimension!\n"));

    for(k=j=0; j<n1; j++)
    	for(i=0;i<j+1; i++,k++)
    	   A[i*n1+j] = A[j*n1+i] = AP[k];
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
	 error(_("No trend order specified!\n"));

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
	      error(_("Trend matrix is singular!\n"));
	for(int i = 0; i < n; i++, pF++) *pF = 1;
	if(trend > 0) {
		//linear terms
	   if(n < d+2)
	     error(_("Trend matrix is singular!\n"));

	   for(int j = 0; j < d; j++) //columns of x
		   for(int k = 0; k < n; k++,pF++) // rows of F
			    *pF = x[j*n+k];
	}
	if(trend > 1) {
	 //quadratic terms
		if( n < (d+1)*(d+2)/2+1 )
	      error(_("Trend matrix is singular!\n"));

		for (int k = 0; k < d; k++)
		  for (int l = k; l < d; l++)
		    for ( int i = 0; i < n; i++, pF++)
			  *pF = x[k*n+i]*x[l*n+i];

	} else if(trend > 2){
	    ERR("Invalid trend order number. Choose either k=1,2.");
	}

}
