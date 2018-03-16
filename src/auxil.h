/**
 * @file        auxil.h
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       R (Interface) Functions
 *
 *
 */
#ifndef AUX_H_
#define AUX_H_

template<class T>
inline T SQR(const T a) {return a*a;}

template<class T>
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

typedef void (*Fp_merger)(double *, double*, int*);
typedef void (*fnCompare_wrap) (double *x, void *data, double *f, int &err);

typedef struct qr_data_s {
	  double *qr, *tau;
	  int  rank, nrow, ncol, ntau;
} qr_data_t, *qr_data;

qr_data qr_dgeqrf(double *a, int m, int n, int &err);

typedef enum  { Chol = 1, SVD = 2, Bunch = 3 } inversion_type;

#ifdef  __cplusplus
extern "C" {
#endif

int  chol2var(double *x, double *z, int nx, double *y );

void fdJacobian(double *x, int dx, double *fval, int m, double *jac, double *fdwork,
	  fnCompare_wrap func, void* data, double eps, int to_negative, int &info);

void dist_X1_X2( double *x, int nrx, int xdim, double *y, int nry, double *d );

void isPositiveDefinite(double *C, int *dim, int *res);
void invMatrix(double *A, int nA, double *res, int &err, inversion_type type = Chol);
void trendfunc(double *x, int lx, int dx, double *f, int model);
//
void Fmatrix(double *x,double *F, int n, int d, int trend);
void triangMat_U_back(double *A, int nA, double *AP, int &info);
void triangMat_U(double *A, int nA, double *AP, int &info);
void set2diag(double *x, double *y, int nx, double s);
int  add2diag(double *x, int nx, double *s);
void splitDiag(double *x, int nx, double *d);
//
void matmult(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z, int &info);
void matmult_trans(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z, int &info);
void matmult_diag(double *x, int nrx, int ncx, double *y, int &info);
void matmult_diag_sqrt(double *x, int nrx, int ncx, double *y, int &info);
void mat_trans(double *y, int ldy, double *x, int ldx, int nrow, int ncol, int &info);
void matrix_col_sum(double *x, int *nx, int *y);

void solveLU(double *A, int nA, double *B, int nB, int &err);
void solveQR(double *X, int nrx, int ncx, double *y, int ncy, int &err);
//
void solveCH(double *X, int n, double *y, int ncoly, int &err);
void solveDSPTRS( double *A, int n, double *B, int nrhs, int &err);
void gsiSolve(double *A, int n, double *B, int nrhs, double *Awork, int &err, inversion_type type);
//
void factorize_chol_L(double *A, int *nA, int &err);
void solve_chol_factorized(double *A, int *nA, double *B, int *nB, int &err);
void solve_chol_triangular(double *A, int *nA, double *B, int *nB, int &err);
//
void Imat(double *x, int n);
//
double norm_x_2( double *x, int nrx, int ncx, int ii);
double norm_2( double *x1, double *x2, int dx1, int dx2, int xdim);
double norm2( double *x1, int n1, double *x2, int n2, int d, int i1, int i2);
//
double denorm(double *x, int n);
double innerProduct(double *x, double *y, int size);
//
void nullSpaceMat(qr_data qr, double *x,int &err);
void qrFree(qr_data qr);
//
void basis0(double *f);
void basis1(double *x, int lx, int dx, double *f);
void basis2(double *x, int lx, int dx, double *f);



#ifdef  __cplusplus
}
#endif

#endif /* AUX_H_ */
