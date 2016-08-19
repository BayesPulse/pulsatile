//
// Test function for getting dataframes from R to C
//


//#define R_NO_REMAP 
//#include <R.h>
//#include <Rinternals.h>
// For boolean variables, consider: #include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "r_interface.h"
#include "pulse_node.h"
#include "birth_death.h"
#include "mcmc.h"

double fitstart; // First time a pulse can occur (10 min. increments).
double fitend;   // Last time a pulse can occur (10 min. increments).
double mmm = 3;


//
// read in arguments and data and convert to types for original decon code
//
SEXP r_interface(SEXP indata,
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

  int i;
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
  double propsd[7];   // Array of proposal variances for MH algorithms
  propsd[4] = Rf_asReal(pv_mean_pulse_mass); 
  propsd[5] = Rf_asReal(pv_mean_pulse_width);
  propsd[2] = Rf_asReal(pv_pulse_mass); 
  propsd[3] = Rf_asReal(pv_pulse_width);
  propsd[6] = Rf_asReal(pv_pulse_location); 
  propsd[0] = Rf_asReal(pv_baseline); 
  propsd[1] = Rf_asReal(pv_halflife);

  Rprintf("proposal variance for location is %f\n", propsd[6]);



  // Set up pulse linklist ------------------
  //   NOTE: Possibly add option to pulse_spec for setting starting values/num
  //   pulses
  double time[9];       // Starting values: pulse locations
  double mass[9];       // Starting values: pulse masses
  double width[9];      // Starting values: pulse widths
  Node_type *list;      
  Node_type *new_node;  
  list = initialize_node();
  mass[0]  = 12.32657558;
  mass[1]  = 8.335093817;
  mass[2]  = 6.169354435;
  mass[3]  = 17.80691618;
  mass[4]  = 17.81203064;
  mass[5]  = 4.180821883;
  mass[6]  = 9.385968508;
  mass[7]  = 3.539176735;
  mass[8]  = 8.860152668;
  width[0] = 26.02652032;
  width[1] = 37.383667;
  width[2] = 37.90962662;
  width[3] = 31.45537546;
  width[4] = 52.64992824;
  width[5] = 42.00314896;
  width[6] = 11.7825812;
  width[7] = 21.11221475;
  width[8] = 100.6173752;
  time[0]  = 139.2555248;
  time[1]  = 293.9643442;
  time[2]  = 400.3647507;
  time[3]  = 654.7396597;
  time[4]  = 788.742805;
  time[5]  = 918.9406528;
  time[6]  = 1156.007127;
  time[7]  = 1335.413248;
  time[8]  = 1408.372033;

  // Initialize nodes and insert into linkedlist
  for (i = 0; i < 9; i++) {

    new_node               = initialize_node();
    new_node->time         = time[i];
    new_node->theta[0]     = mass[i];
    new_node->theta[1]     = width[i];
    new_node->mean_contrib = (double *)calloc(nobs, sizeof(double));
    mean_contribution(new_node, ts, parms, nobs);
    insert_node(new_node, list);

  }

  //print_list(list);



  // R lists for chains ------------------
  // Common parameters -- 1 record per saved iteration
  SEXP common1 = Rf_protect(Rf_allocMatrix(REALSXP, iters/nthin, 8));
  // Pulse-specific parameters -- 1 record per pulse per iteration
  //   list object defined here, each iteration will differ in lenght, so the
  //   entries are defined in linklistv3
  SEXP parm1   = Rf_protect(Rf_allocVector(VECSXP, iters/nthin));


  GetRNGstate();
  mcmc(list, parms, ts, iters, nobs, nthin, priors, common1, parm1, propsd);
  PutRNGstate();



  // Combine chains for output
  SEXP chains = Rf_protect(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(chains, 0, common1);
  SET_VECTOR_ELT(chains, 1, parm1);




  // Free memory
  deallocate_data(ts, nobs);
  destroy_list(list);
  free(priors->fe_mean);
  free(priors->fe_variance);
  free(priors->re_sdmax);
  free(priors);
  free(parms->re_sd);
  free(parms);
  Rf_unprotect(3);

  //return Rf_ScalarReal(nthin);
  //return Rf_asVector(priors);
  return chains;

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

