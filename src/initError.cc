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

int maxErrors = MAX_ERRORS, lastError = 0;
int maxWarnings = MAX_WARNINGS, lastWarning = 0;

char R_MSG_BUFFER[1000]="";
char C_MSG_BUFFER[1000]="";
char ERROR_STRING[maxErrorStr]="";
char ERROR_LOC[nErrorLocations]="";


/* error and warnings containers*/
error_queue errContainer;
warning_queue wrrContainer;

/* get error message from id*/
void errorMSG(int err, char* m) {
  if (PL > 0) {
      PRINTF("\n\n #### C ERROR MSG block #### \n");
      PRINTF("error code %d\n", err);
  }

  if (err >= ERROR_NUM_MIN && err <= ERROR_NUM_MAX) err = ERROR_NUM_MIN;

  switch (err) {
      case NO_ERROR : return;
      case NaN_ERROR: std::strcpy(m,"NaNs produced"); break;
      case LAPACK_ERROR: std::strcpy(m,"Lapack function returned error code"); break;
      case LAPACK_QR_ERROR: std::strcpy(m,"Lapack function error in routine qr_dgeqrf!"); break;
      case LAPACK_PMAT_ERROR: std::strcpy(m,"Error in routine nullSpaceMat!"); break;
      case LAPACK_FACTORIZE_ERROR:  std::strcpy(m,"Error in routine dpotrf!"); break;
      case LAPACK_SOLVE_ERROR:  std::strcpy(m,"Error in routine dpotrs!"); break;
      case LAPACK_INVERSION_ERROR:  std::strcpy(m,"Error in routine invMatrix!"); break;
      default:
    	 std::strcpy(m,"default error: ");
         break;
  }
}

void warningMSG(int wrr, char* m) {
  switch(wrr) {
  case NO_WARNING: return;
  case GENERIC_WARNING: std::strcpy(m,"Generic warning"); break;
  default:
	  std::strcpy(m,"default warning: ");
    break;
  }
}

void printMSG(int err, const char *msg , int line, const char *file){
   char MSG[maxErrorStr];
   errorMSG(err, ERROR_STRING); /* errors and warnings together */
   if(err < WARNING_NUM_MIN)
     std::sprintf(ERROR_LOC, "Error in file %s at line %d. \n ", file, line);
   else
	 std::sprintf(ERROR_LOC, "Warning in file %s at line %d. \n ", file, line);
   std::sprintf(MSG, "%s\n  ... %s\n ... %s\n", ERROR_LOC, ERROR_STRING, msg);
   PRINT(" %s ", MSG);
}

void incError() {
   ++lastError;
   if(lastError==MAX_ERRORS) {
	   std::strcpy(ERROR_LOC,__FILE__);
      cleanErrors();
      ERR(" *** Error buffer overflow! *** ");
   } else if(lastError>MAX_ERRORS/2) {
	   std::strcpy(ERROR_LOC,__FILE__);
     WRR(" *** Risk of error buffer overflow: too many errors/warnings ***.\n");
   }
}

void incWarning() {
   ++lastWarning;
   if(lastWarning==MAX_WARNINGS) {
	   std::strcpy(ERROR_LOC,__FILE__);
      cleanWarnings();
      ERR(" *** Warning buffer overflow! *** ");
   } else if(lastWarning>MAX_WARNINGS/2) {
	   std::strcpy(ERROR_LOC,__FILE__);
     WRR(" *** Risk of warning buffer overflow: too many errors/warnings ***.\n");
   }
}

void setError(int id, const char* msg, int line, const char* file, int code )
{
    errContainer[lastError] = Calloc(1, errorStatus_s);
    errContainer[lastError]->id = id;
    errContainer[lastError]->line = line;
    errContainer[lastError]->isError = 1;
    errContainer[lastError]->numErr = lastError;
    errContainer[lastError]->code = code;

    strcpy(errContainer[lastError]->file,file);
    strcpy(errContainer[lastError]->msg,msg);

    incError();

}

void setWarning(int id, const char* msg, int line, const char* file, int code)
{
    wrrContainer[lastWarning] = Calloc(1, warningStatus_s);
    wrrContainer[lastWarning]->id = id;
    wrrContainer[lastWarning]->line = line;
    wrrContainer[lastWarning]->isWarning = 1;
    wrrContainer[lastWarning]->numWrr = lastWarning;
    wrrContainer[lastWarning]->code = code;

    std::strcpy(wrrContainer[lastWarning]->file,file);
    std::strcpy(wrrContainer[lastWarning]->msg,msg);
    incWarning();

}

void cleanErrors() {
  if(!lastError)
     return;
   while(lastError > 0) {
	Free(errContainer[lastError-1]);
	--lastError;
   }
}

void cleanWarnings() {
  if(!lastWarning)
    return;
  while(lastWarning > 0) {
        Free(wrrContainer[lastWarning-1]);
        --lastWarning;
   }
}

SEXP emptyErrorCache() {
  SEXP R_msg,R_list,R_cls;

  char MSG[maxErrorStr];
  errorStatus *err = NULL;

  PROTECT(R_list = allocVector(VECSXP, lastError));
  PROTECT(R_cls = allocVector(STRSXP, 1));

  SET_STRING_ELT(R_cls, 0, mkChar("error"));
  setAttrib(R_list, R_ClassSymbol, R_cls);

  for(int j=0; j < lastError; j++) {
    err = errContainer[j];
    std::sprintf(MSG,"Error %d, code %d: %s in file `%s' at line %d",err->numErr,err->id,err->msg,err->file,err->line);
    PROTECT(R_msg = allocVector(STRSXP, 1));
    SET_STRING_ELT(R_msg, 0, mkChar(MSG));
    SET_VECTOR_ELT(R_list, j, R_msg);
    UNPROTECT(1);
  }

  if( (lastError>0) )
     cleanErrors();

  UNPROTECT(2);
  return R_list;
}


SEXP emptyWarningCache() {
  SEXP R_msg,R_list,R_cls;

  char MSG[maxErrorStr];
  warningStatus *wrr = NULL;

  PROTECT(R_list = allocVector(VECSXP, lastWarning));
  PROTECT(R_cls = allocVector(STRSXP, 1));

  SET_STRING_ELT(R_cls, 0, mkChar("warning"));
  setAttrib(R_list, R_ClassSymbol, R_cls);

  for(int j=0; j < lastWarning; j++) {
          wrr = wrrContainer[j];
          std::sprintf(MSG,"Warning %d, code %d: %s in file `%s' at line %d",wrr->numWrr,wrr->id,wrr->msg,wrr->file,wrr->line);
          PROTECT(R_msg = allocVector(STRSXP, 1));
          SET_STRING_ELT(R_msg, 0, mkChar(MSG));
          SET_VECTOR_ELT(R_list, j, R_msg);
          UNPROTECT(1);
  }

  if(lastWarning > 0)
     cleanWarnings();

  UNPROTECT(2);
  return R_list;
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



