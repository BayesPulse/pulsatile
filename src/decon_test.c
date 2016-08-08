//
// Test function for getting dataframes from R to C
//


#define R_NO_REMAP 
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "decon_test.h"



//
// read in arguments and data and convert to types for original decon code
//
void decon_input(SEXP indata,
                 SEXP model,
                 SEXP thin,
                 SEXP iterations
                 ) {

  double **ts; // Pointer to Nx2 matrix containing data
  int nobs; // n observations in data (*N)

  // Number of obs in series
  nobs = Rf_xlength(getListElement(indata, "concentration"));

  // Data SEXP to C array
  ts = convert_data(indata, nobs);

  // Algo arguments SEXP to C;
  long iters  = Rf_asInteger(iterations); 
  int nthin = Rf_asInteger(thin);


  //Print to check values
  Rprintf("There are %ld observations in the dataset\n", nobs);
  Rprintf("This data is of type (model arg not yet implemented)");
  Rprintf("The algo was run for %ld iterations\n", iters);
  Rprintf("Every %ld th sample was saved in the output chain", nthin);
  // Read the time series into memory                          
  Rprintf("And the dataset looks like this:\n");
  for (int i = 0; i < nobs; i++) {
    Rprintf("time = %lf, conc = %lf\n", ts[i][0], ts[i][1]);
  }


  // Free memory
  deallocate_data(ts, nobs);

  return ;

}




//
// convert data from SEXP to C array
//
double **convert_data(SEXP indata, int nrow) { 

  // Declare variables
  int i = 0;
  double **data; // Nx2 matrix of time and log-concentration

  // Allocate memory for matrix of data
  data = (double **)calloc(nrow, sizeof(double *));
  for (i = 0; i < nrow; i++) {
    data[i] = (double *)calloc(2, sizeof(double));
  }

  // Read the time series into memory                          
  for (i = 0; i < nrow; i++) {
    data[i][0] = INTEGER(getListElement(indata, "time"))[i];
    data[i][1] = REAL(getListElement(indata, "concentration"))[i];
    data[i][1] = log(data[i][1]);
  }

	return data;
	
}


void deallocate_data(double **data, int nrow) {

  int i = 0;

  for (i = 0; i < nrow; i++) {
    free(data[i]);
  }

  free(data);
}



//
// Fn for accessing list elements
//
SEXP getListElement(SEXP list, const char *str) {

  SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);

  for (int i = 0; i < Rf_xlength(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }

  return elmt;

}



//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

