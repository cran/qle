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

int PL = 1;

char R_MSG_BUFFER[100]="";
char C_MSG_BUFFER[100]="";
char ERROR_LOC[nErrorLocations]="";

/* get error message from id*/
void errorMSG(int errId, char* m) {
  if (PL > 0) {
      PRINTF("\n\n #### C ERROR MSG block #### \n");
      PRINTF("error code %d\n", errId);
  }

  switch (errId) {
      case NO_ERROR : return;
      case NaN_ERROR: std::strcpy(m,"NaNs produced"); break;
      case LAPACK_ERROR: std::strcpy(m,"Lapack function returned error code"); break;
      case LAPACK_QR_ERROR: std::strcpy(m,"Lapack function error in routine qr_dgeqrf!"); break;
      case LAPACK_PMAT_ERROR: std::strcpy(m,"Error in routine nullSpaceMat!"); break;
      case LAPACK_FACTORIZE_ERROR:  std::strcpy(m,"Error in routine dpotrf!"); break;
      case LAPACK_SOLVE_ERROR:  std::strcpy(m,"Error in routine dpotrs!"); break;
      case LAPACK_INVERSION_ERROR:  std::strcpy(m,"Error in routine invMatrix!"); break;
      default: std::strcpy(m,"unknown error code: "); break;
  }
}

void warningMSG(int wrr, char* m) {
  switch(wrr) {
   case NO_WARNING: return;
   default: std::strcpy(m,"generic warning code: ");  break;
  }
}

void printMSG(int err, const char *msg , int line, const char *file){
   char MSG[maxErrorStr];
   errorMSG(err, C_MSG_BUFFER); /* errors and warnings together */

   if(err < WARNING_NUM_MIN)
     std::sprintf(ERROR_LOC, "Error in file %s at line %d. \n ", file, line);
   else
	 std::sprintf(ERROR_LOC, "Warning in file %s at line %d. \n ", file, line);

   std::sprintf(MSG, "%s\n  ... %s\n ... %s\n", ERROR_LOC, C_MSG_BUFFER, msg);
   PRINT(" %s ", MSG);
}

void printMatrix(const char ch[], const double *mat, int *irows, int *icols) {
	int j, k, rows, cols;

	rows = *irows;
	cols = *icols;
        PRINT("%s: (%d x %d)\n", ch, rows, cols);

	for(j = 0; j < rows; j++) {
		for(k = 0; k < cols; k++) {
            PRINT("\t%3.16f", mat[k*rows+j]);
		}
        PRINT("\n");
	}
    PRINT("\n");
}

void printVector(const char ch[], const double *vec, int *lx) {
	int j;
	PRINT("%s: length=[ %d ]\n", ch, *lx);
	for(j = 0; j < *lx; j++) {
		PRINT("\t%3.16f", vec[j]);
		if( (j>0 ) && (j % 10 == 0)) Rprintf("\n");
	}
	PRINT("\n");
}

void print_R_matrix( SEXP A, const std::string& matrix_name) {
     std::ostringstream oss;
     if(!isMatrix(A))
    	 error("print_matrix failed: Matrix type expected!\n");
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
     PRINT("%s\n",oss.str().c_str());
}


void print_R_vector( SEXP v, const std::string& vector_name) {
     std::ostringstream oss;
     char numstr[21];
     if(!isVector(v))
    	 error("print_vector failed: Vector expected!\n");
     std::size_t lv = Rf_length(v);
     oss << "Vector: " << vector_name << "\n";
     for(std::size_t i=0;i<lv; i++ ) {
    	 std::sprintf(numstr, "%3.12f", REAL(v)[i]);
         oss << numstr << " ";
     }
     PRINT("%s\n",oss.str().c_str());
 }



