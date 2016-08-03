//
// Test function for getting dataframes from R to C
//


#define R_NO_REMAP 
#include <R.h>
#include <Rinternals.h>
#include "test_inout.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>


// Fn for testing reading in data
SEXP decon(SEXP indata,
           SEXP model,
           SEXP iterations,
           SEXP prior_pulse_mass_mean,
           SEXP prior_pulse_mass_var,
           SEXP prior_pulse_width_mean,
           SEXP prior_pulse_width_var,
           SEXP prior_pulse_location_gamma,
           SEXP prior_pulse_location_range,
           SEXP prior_pulse_location_count,
           SEXP prior_max_sd_mass,
           SEXP prior_max_sd_width,
           SEXP prior_baseline_mean,
           SEXP prior_baseline_var,
           SEXP prior_halflife_mean,
           SEXP prior_halflife_var,
           SEXP prior_error_alpha,
           SEXP prior_error_beta,
           SEXP sv_pulse_mass_mean,
           SEXP sv_pulse_mass_sd,
           SEXP sv_pulse_width_mean,
           SEXP sv_pulse_width_sd,
           SEXP sv_baseline_mean,
           SEXP sv_halflife_mean,
           SEXP sv_error_var,
           SEXP pv_mean_pulse_mass,
           SEXP pv_mean_pulse_width,
           SEXP pv_pulse_mass,
           SEXP pv_pulse_width,
           SEXP pv_pulse_location,
           SEXP pv_baseline,
           SEXP pv_halflife) {

  double **data;
  data = convert_data(indata);
  int i = 0;
  int nrow = Rf_xlength(getListElement(indata, "concentration"));

  for (i = 0; i < nrow; i++) {
    Rprintf("time = %lf, conc = %lf\n", data[i][0], data[i][1]);
  }
  
  for (i = 0; i < nrow; i++) {
    free(data[i]);
  }
  free(data);

  return(Rf_ScalarReal(1));

}


// Fn for testing reading in data
double **convert_data(SEXP indata) { 

  // Declare variables
  int nrow;
  int i = 0;
  double **data; // Nx2 matrix of time and log-concentration

  // Number of obs in series
  nrow = Rf_xlength(getListElement(indata, "concentration"));

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

  //free(nrow);
  //free(i);

	return data;
	
}

// Fn for testing reading in data
SEXP testc(SEXP indata) { //, SEXP specification) {

  // Declare variables
  int nrow;
  int i = 0;
  double **data;   // Nx2 matrix of time and log-concentration

  // Protect R data file
	//indata = Rf_protect(Rf_coerceVector(indata, VECSXP));
	//inspec = PROTECT(coerceVector(specification, VECSXP));

  // Number of obs in series
  nrow = Rf_xlength(getListElement(indata, "concentration"));

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
  }

  // Unprotect R objects and free memory
	//Rf_unprotect(1);
  // Free data memory 
  for (i = 0; i < nrow; i++) {
    free(data[i]);
  }
  free(data);

	return Rf_ScalarReal(nrow);
	
}

// Fn for testing reading in arguments
SEXP testspec(SEXP inspec) { 

  // Declare variables
  int n_firstlvl;
  //int i = 0;
  //double **speclist;   // Nx2 matrix of time and log-concentration

  // Protect R data file
	//inspec = Rf_protect(Rf_coerceVector(inspec, VECSXP));
	//inspec = Rf_protect(coerceVector(specification, VECSXP));

  // Number of obs in series
  n_firstlvl = Rf_xlength(inspec); //getListElement(speclist, "concentration"));

  // Allocate memory for matrix of data
  //data = (double **)calloc(nrow, sizeof(double *));
  //for (i = 0; i < nrow; i++) {
  //  data[i] = (double *)calloc(2, sizeof(double));
  //}

  //// Read the time series into memory                          
  //for (i = 0; i < nrow; i++) {
  //  data[i][0] = INTEGER(getListElement(indata, "time"))[i];
  //  data[i][1] = REAL(getListElement(indata, "concentration"))[i];
  //  data[i][1] = log(data[i][1]);
  //  Rprintf("time = %lf, conc = %lf\n", data[i][0], data[i][1]);
  //}

  // Unprotect R objects and free memory
	//Rf_unprotect(1);
  // Free data memory 
  //for (i = 0; i < nrow; i++) {
  //  free(speclist[i]);
  //}
  //free(speclist);

	return Rf_ScalarReal(n_firstlvl);
	
}

// Fn for showing arg types in any data type
SEXP showArgs1(SEXP largs) {

  int i, nargs = Rf_xlength(largs);
  Rcomplex cpl;
  SEXP el, names = Rf_getAttrib(largs, R_NamesSymbol);
  const char *name;

  for(i = 0; i < nargs; i++) {

    el = VECTOR_ELT(largs, i);
    name = Rf_isNull(names) ? "" : CHAR(STRING_ELT(names, i));

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

