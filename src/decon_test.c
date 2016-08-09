//
// Test function for getting dataframes from R to C
//


#define R_NO_REMAP 
#include <R.h>
#include <Rinternals.h>
//#include <stdlib.h>
//#include <math.h>
//#include <stdio.h>
//#include <time.h>
#include "decon_test.h"

double fitstart; // First time a pulse can occur (10 min. increments).
double fitend;   // Last time a pulse can occur (10 min. increments).


//
// read in arguments and data and convert to types for original decon code
//
SEXP decon_input(SEXP indata,
                 SEXP model,
                 SEXP thin,
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
                 SEXP prior_error_beta
                 ) {

  double **ts; // Pointer to Nx2 matrix containing data
  int nobs; // n observations in data (*N)

  // Number of obs in series
  nobs = Rf_xlength(getListElement(indata, "concentration"));

  // Data SEXP to C array
  ts = convert_data(indata, nobs);

  // Algo arguments SEXP to C;
  long iters  = Rf_asInteger(iterations); 
  long nthin = Rf_asInteger(thin);


  //Print to check values
  Rprintf("There are %d observations in the dataset.\n", nobs);
  Rprintf("This data is of type (model arg not yet implemented).\n");
  Rprintf("The algo was run for %ld iterations.\n", iters);
  Rprintf("Every %ldth sample was saved in the output chain.\n", nthin);
  // Read the time series into memory                          
  Rprintf("And the dataset looks like this:\n");
  for (int i = 0; i < nobs; i++) {
    Rprintf("time = %lf, conc = %lf\n", ts[i][0], ts[i][1]);
  }


  fitend   = ts[nobs - 1][0] + ts[0][0] * 2; // search 2 units farther in time
  fitstart = -ts[0][0] * 4;                // search 4 units in the past

  // Set up priors structure ------------------
  Priors *priors;                 // Prior parameters data structure
  priors                 = (Priors *)calloc(1, sizeof(Priors));
  priors->re_sdmax       = (double *)calloc(2, sizeof(double));
  priors->fe_variance    = (double *)calloc(2, sizeof(double));
  priors->fe_mean        = (double *)calloc(2, sizeof(double));
  priors->meanbh[0]      = Rf_asReal(prior_baseline_mean);    // priormub;
  priors->meanbh[1]      = Rf_asReal(prior_halflife_mean);    // priormuh;
  priors->varbh[0]       = Rf_asReal(prior_baseline_var);     // priorvarb;
  priors->varbh[1]       = Rf_asReal(prior_halflife_var);     // priorvarh;
  priors->fe_mean[0]     = Rf_asReal(prior_pulse_mass_mean);  // priormu1;
  priors->fe_mean[1]     = Rf_asReal(prior_pulse_width_mean); // priormu2;
  priors->fe_variance[0] = Rf_asReal(prior_pulse_mass_var);   // priorvar1;
  priors->fe_variance[1] = Rf_asReal(prior_pulse_width_var);  // priorvar2;
  priors->re_sdmax[0]    = Rf_asReal(prior_max_sd_mass);      // priora1;
  priors->re_sdmax[1]    = Rf_asReal(prior_max_sd_width);     // priora2;
  priors->err_alpha      = Rf_asReal(prior_error_alpha);      // prioralpha;
  priors->err_beta       = Rf_asReal(prior_error_beta);       // priorbeta;
  priors->gamma          = Rf_asReal(prior_pulse_location_gamma);   // priorgamma;
  priors->range          = Rf_asReal(prior_pulse_location_range);   // priorrange;

  // Free memory
  deallocate_data(ts, nobs);

  return Rf_ScalarReal(nthin);
  //return Rf_asVector(priors);

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


//
// Deallocate data set array
//
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

  SEXP elmt  = R_NilValue; 
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);

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

