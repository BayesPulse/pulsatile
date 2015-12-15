//
// Test function for getting dataframes from R to C
//


#include "test_inout.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>


// Fn for testing reading in data
SEXP testc(SEXP indata) {

  // Declare variables
	//SEXP ab;
	//SEXP mean_conc;
	//double mean_conc_c = 0;
	//int xnrow = 0;
  //double *conc_dat;
  int nrow;
  int i = 0;
  double **data;   // Nx2 matrix of time and log-concentration

	// Think this is how should pass list objects back to R 
	//common = PROTECT(allocVector(VECSXP, common_c));
	//pulse  = PROTECT(allocVector(VECSXP, pulse_c));
	//nrow = PROTECT(coerceVector(nrow, INTSXP));
	//mean_conc = PROTECT(allocVector(REALSXP, mean_conc_c));
	//conc_dat = PROTECT(coerceVector(data, VECSXP));
  //xnrow = REAL(nrow);
  
  // Protect R data file
	indata = PROTECT(coerceVector(indata, VECSXP));

  // Number of obs in series
  nrow = LENGTH(getListElement(indata, "concentration"));

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
    Rprintf("time = %lf, conc = %lf\n", data[i][0], data[i][1]);
    //Rprintf("conc = %lf\n", REAL(getListElement(indata, "concentration"))[i]);
  }
	//get 'length' of data (whatever length is for this object)
	//nrow = LENGTH(data);
  //Rprintf("Data is %d\n", TYPEOF(data));
  //nrow = LENGTH(getListElement(data, "concentration"));
  //Rprintf("Len of conc var is  %d\n", nrow);


  //for (int i = 0; i < nrow; i++) {
  //  Rprintf("conc_dat = %d\n", conc_dat[i]);
  //  //conc_dat[i] = asReal(data[i]);
  //}

	UNPROTECT(1);
  // Free data memory 
  for (i = 0; i < nrow; i++) {
    free(data[i]);
  }
  free(data);
  //free(nrow);

	return(ScalarReal(nrow));
	

	//a = PROTECT(coerceVector(a, REALSXP));
	//b = PROTECT(coerceVector(b, REALSXP));
	//na = length(a); 
	//nb = length(b); 
	//nab = na + nb - 1;
	//ab = PROTECT(allocVector(REALSXP, nab));
	//xa = REAL(a); 
	//xb = REAL(b); 
	//xab = REAL(ab);
	//for(int i = 0; i < nab; i++) {
	//	xab[i] = 0.0;
	//}
	//for(int i = 0; i < na; i++) {
	//	for(int j = 0; j < nb; j++) { 
	//		xab[i + j] += xa[i] * xb[j];
	//	}
	//}
	//UNPROTECT(3);
	//return ab;

}


// Fn for showing arg types in any data type
SEXP showArgs1(SEXP largs) {

    int i, nargs = LENGTH(largs);
    Rcomplex cpl;
    SEXP el, names = getAttrib(largs, R_NamesSymbol);
    const char *name;

    for(i = 0; i < nargs; i++) {
  el = VECTOR_ELT(largs, i);
  name = isNull(names) ? "" : CHAR(STRING_ELT(names, i));
  switch(TYPEOF(el)) {
  case REALSXP:
      Rprintf("[%d] '%s' %f\n", i+1, name, REAL(el)[0]);
      break;
  case LGLSXP:
  case INTSXP:
      Rprintf("[%d] '%s' %d\n", i+1, name, INTEGER(el)[0]);
      break;
  case CPLXSXP:
      cpl = COMPLEX(el)[0];
      Rprintf("[%d] '%s' %f + %fi\n", i+1, name, cpl.r, cpl.i);
      break;
  case STRSXP:
      Rprintf("[%d] '%s' %s\n", i+1, name,
        CHAR(STRING_ELT(el, 0)));
      break;
  default:
      Rprintf("[%d] '%s' R type\n", i+1, name);
  }
    }
    return R_NilValue;
}


// Fn for accessing list elements
SEXP getListElement(SEXP list, const char *str) {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    for (int i = 0; i < length(list); i++)
  if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
     elmt = VECTOR_ELT(list, i);
     break;
  }
    return elmt;
}



//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

