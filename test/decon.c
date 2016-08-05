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


// includes/defs from deconvolution_main.c
#include "format_data.h"
#include "hash.h"
#include "linklistv3.h"
#include "birthdeath_strauss.h"
#define EPS 1.0e-12 // error for testing equality with 0

// Global variables from deconvolution_main.c
double M; // Precision used in RNGs. Set to 32 or 64 based on operating system.
int mmm = 3;     // Order statistic to use for order stat prior
double fitstart; // First time a pulse can occur (10 min. increments).
double fitend;   // Last time a pulse can occur (10 min. increments).

//
// Primary analysis function 
//
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

  // Print to verify data read in
  //for (i = 0; i < nrow; i++) {
  //  Rprintf("time = %lf, conc = %lf\n", data[i][0], data[i][1]);
  //}
  






  for (i = 0; i < nrow; i++) {
    free(data[i]);
  }
  free(data);

  return(Rf_ScalarReal(1));

}


//
// Fn for testing reading in data
//
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

