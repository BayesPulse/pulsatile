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
                 SEXP pv_halflife
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
  Rprintf("And the first 10 obs of the dataset looks like this:\n");
  for (int i = 0; i < 10; i++) {
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

  Rprintf("The fit start is %f and the fit end is at time %f\n", fitstart, fitend);
  Rprintf("The prior on the baseline has mean %f and variance %f\n", priors->meanbh[0], priors->varbh[0]);

  // Set up parms structure ------------------
  Common_parms *parms;            // Common parameters data structure
  parms           = (Common_parms *)calloc(1, sizeof(Common_parms));
  parms->re_sd    = (double *)calloc(2, sizeof(double));
  parms->nprior   = Rf_asReal(prior_pulse_location_count); // priorr;
  parms->theta[0] = Rf_asReal(sv_pulse_mass_mean);   // svmu1;
  parms->theta[1] = Rf_asReal(sv_pulse_width_mean);  // svmu2;
  parms->md[0]    = Rf_asReal(sv_baseline_mean);     // svbase;
  parms->md[1]    = Rf_asReal(sv_halflife_mean);     // svhalf;
  parms->sigma    = Rf_asReal(sv_error_var);         // svevar; // error variance
  parms->lsigma   = log(parms->sigma);    // log of error variance
  parms->re_sd[0] = Rf_asReal(sv_pulse_mass_sd);  // svsig1; 
  parms->re_sd[1] = Rf_asReal(sv_pulse_width_sd); // svsig2; 
  // note: other re_sd parms (max and pv) are entered and used on the natural
  // scale, so for consistency, we enter these on the same scale.  Functions
  // other than draw_re_sd, use re_sd[j] on the log scale so we log xform them
  // here.                                                                   
  parms->decay    = log(2) / parms->md[1]; // calculate decay rate from halflife


  // Set up proposal variances ------------------
  double propvar[7];   // Array of proposal variances for MH algorithms
  propvar[4] = Rf_asReal(pv_mean_pulse_mass); 
  propvar[5] = Rf_asReal(pv_mean_pulse_width);
  propvar[2] = Rf_asReal(pv_pulse_mass); 
  propvar[3] = Rf_asReal(pv_pulse_width);
  propvar[6] = Rf_asReal(pv_pulse_location); 
  propvar[0] = Rf_asReal(pv_baseline); 
  propvar[1] = Rf_asReal(pv_halflife);

  Rprintf("proposal variance for location is %f\n", propvar[6]);

  // Free memory
  deallocate_data(ts, nobs);
  free(priors->fe_mean);
  free(priors->fe_variance);
  free(priors->re_sdmax);
  free(priors);
  free(parms->re_sd);
  free(parms);

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

