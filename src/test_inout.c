//
// Test function for getting dataframes from R to C
//
//
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

SEXP test_inout_(SEXP data) {

	//int na, nb, nab;
	//double *xa, *xb, *xab;
	//SEXP ab;

	//SEXP mean_conc;
	//double mean_conc_c = 0;
	//int xnrow = 0;
  int nrow;

	// Think this is how should pass list objects back to R 
	//common = PROTECT(allocVector(VECSXP, common_c));
	//pulse  = PROTECT(allocVector(VECSXP, pulse_c));
	data = PROTECT(coerceVector(data, VECSXP));
	//nrow = PROTECT(coerceVector(nrow, INTSXP));
	//mean_conc = PROTECT(allocVector(REALSXP, mean_conc_c));
  
  //xnrow = REAL(nrow);

	//get 'length' of data (whatever length is for this object)
	nrow = TRUELENGTH(data);

	UNPROTECT(1);
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
