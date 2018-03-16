/**
 * @file        initError.cc
 * @author      M. Baaske
 * @date        11/03/2013
 * @brief       Initialize error functions
 *
 * Explanation:
 *
 */

#include "error.h"

/*
 * PL = 1   : show warnings
 * PL <= 10 : do not show warnings from `gsiSolve` and `gsiInv` except NA warnings
 */
int PL = 1;

char C_MSG_BUFFER[1000] = "";

/* get error message from id*/
void errorMSG(int errId, const char *name) {
    char msg[100] = "";
	switch (errId) {
      case NO_ERROR : return;
      case NaN_ERROR: std::strcpy(msg,"NaNs detected."); break;
      case LAPACK_ERROR: std::strcpy(msg,"Generic Lapack function error."); break;
      case LAPACK_QR_ERROR: std::strcpy(msg,"Lapack QR decomposition failed."); break;
      case LAPACK_PMAT_ERROR: std::strcpy(msg,"Error constructing kernel (null-space) matrix."); break;
      case LAPACK_FACTORIZE_ERROR:  std::strcpy(msg,"Lapack error in factorization."); break;
      case LAPACK_SOLVE_ERROR:  std::strcpy(msg,"Lapack error solving linear system of equations."); break;
      case LAPACK_INVERSION_ERROR:  std::strcpy(msg,"Lapack error in inverting matrix."); break;
      default: std::strcpy(msg,"Unknown error code: "); break;
  }
  std::sprintf(C_MSG_BUFFER,"error in function `%s` (code=%d): %s.\n", name, errId, msg);
}

void warningMSG(int wrrId, const char *name) {
	 char msg[100] = "";

	 switch(wrrId) {
	   case NO_WARNING: return;
	   case NaN_WARNING: std::strcpy(msg,"NaNs detected."); break;
	   case LAPACK_WARNING: std::strcpy(msg,"Lapack error in factorization."); break;
	   case POSDEF_WARNING: std::strcpy(msg,"Matrix probably not (semi)positive definite."); break;
	   default: std::strcpy(msg,"generic warning code: ");  break;
	  }
	 std::sprintf(C_MSG_BUFFER,"warning in function `%s` (code=%d): %s.\n", name, wrrId, msg);
}

void printMatrix(const char ch[], const double *mat, int *irows, int *icols) {
	int j=0, k=0, rows=*irows, cols=*icols;
    Rprintf("%s: (%d x %d)\n", ch, rows, cols);

	for(j = 0; j < rows; j++) {
		for(k = 0; k < cols; k++) {
            Rprintf("\t%3.16f", mat[k*rows+j]);
		}
        Rprintf("\n");
	}
    Rprintf("\n");
}

void printVector(const char ch[], const double *vec, int *lx) {
	int j=0;
	Rprintf("%s: length=[ %d ]\n", ch, *lx);
	for(j = 0; j < *lx; j++) {
		Rprintf("\t%3.16f", vec[j]);
		if( (j>0 ) && (j % 10 == 0)) Rprintf("\n");
	}
	Rprintf("\n");
}

void print_R_matrix( SEXP A, const std::string& matrix_name) {
     std::ostringstream oss;
     if(!isMatrix(A))
      Rf_error("print_matrix failed: Matrix type expected!\n");
     int *dimX = GET_DIMS(A);
     std::size_t size1 = dimX[0], size2=dimX[1];
     char numstr[21]; // enough to hold all numbers up to 64-bits
      oss << "Matrix: " << matrix_name << "\n"
          << "[size m,size n] = " << "[" << size1 <<","<< size2 <<"]" << "\n";
	 for (std::size_t i = 0; i < size1; i++){
          oss << "[" << i << "] ";
                for (std::size_t j = 0; j < size2; j++) {
                	std::sprintf(numstr, "\t%3.16f", REAL(A)[i+size1*j]);
                    oss << numstr << " ";
                }
                oss << "\n";
	 }
     oss << "\n";
     Rprintf("%s\n",oss.str().c_str());
}


void print_R_vector( SEXP v, const std::string& vector_name) {
     std::ostringstream oss;
     char numstr[21];
     if(!isVector(v))
       Rf_error("print_vector failed: Vector expected!\n");
     std::size_t lv = Rf_length(v);
     oss << "Vector: " << vector_name << "\n";
     for(std::size_t i=0;i<lv; i++ ) {
    	 std::sprintf(numstr, "%3.12f", REAL(v)[i]);
         oss << numstr << " ";
     }
     Rprintf("%s\n",oss.str().c_str());
 }



