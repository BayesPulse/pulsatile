//
// Primary interface between R and C code
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


//
// read in arguments and data and convert to types for original decon code
//
SEXP decon_r_interface(SEXP indata,
                       SEXP model,
                       SEXP thin,
                       SEXP burnin,
                       SEXP iterations,
                       SEXP inverbose,
                       SEXP strauss_location_prior,
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
                       SEXP pv_indiv_pulse_mass,
                       SEXP pv_indiv_pulse_width,
                       SEXP pv_sd_pulse_mass,
                       SEXP pv_sd_pulse_width,
                       SEXP pv_pulse_location,
                       SEXP pv_baseline,
                       SEXP pv_halflife,
                       SEXP pv_etamass,
                       SEXP pv_etawidth
                       ) {

  //int i;
  double **ts; // Pointer to Nx2 matrix containing data
  int nobs; // n observations in data (*N)

  // Number of obs in series
  nobs = Rf_xlength(getListElement(indata, "concentration"));

  // Data SEXP to C array
  ts = convert_data(indata, nobs);

  // Algo arguments SEXP to C;
  long iters  = Rf_asInteger(iterations); 
  int verbose = Rf_asInteger(inverbose); 
  int nthin   = Rf_asInteger(thin);
  int strauss = Rf_asInteger(strauss_location_prior);
  int nburnin = Rf_asInteger(burnin);
  int out_len = (iters - nburnin) / nthin;
  

  //Print to check values
  //Rprintf("There are %d observations in the dataset.\n", nobs);
  //Rprintf("This data is of type (model arg not yet implemented).\n");
  //Rprintf("The algo was run for %ld iterations.\n", iters);
  //Rprintf("Every %ldth sample was saved in the output chain.\n", nthin);
  //// Read the time series into memory                          
  //Rprintf("And the first 10 obs of the dataset looks like this:\n");
  //for (int i = 0; i < 10; i++) {
  //  Rprintf("time = %lf, conc = %lf\n", ts[i][0], ts[i][1]);
  //}

  fitend   = ts[nobs - 1][0] + ts[0][0] * 2; // search 2 units farther in time
  fitstart = -(ts[1][0] - ts[0][0]) * 4;     // search 4 units in the past

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
  priors->orderstat      = 3;                                 // hard-coded;
  priors->gamma          = Rf_asReal(prior_pulse_location_gamma);   // priorgamma;
  priors->range          = Rf_asReal(prior_pulse_location_range);   // priorrange;

  //Rprintf("The fit start is %f and the fit end is at time %f\n", fitstart, fitend);
  //Rprintf("The prior on the baseline has mean %f and variance %f\n", priors->meanbh[0], priors->varbh[0]);

  // Set up parms structure ------------------
  Common_parms *parms;            // Common parameters data structure
  parms           = (Common_parms *)calloc(1, sizeof(Common_parms));
  //parms->re_sd    = (double *)calloc(2, sizeof(double));
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
  double propsd[11];   // Array of proposal variances for MH algorithms
  propsd[4]  = Rf_asReal(pv_indiv_pulse_mass);  // SD of proposal distr for random effects (individ. pulse masses) -- formerly pv_mean_pulse_mass
  propsd[5]  = Rf_asReal(pv_indiv_pulse_width); // SD of proposal distr for random effects (individ. pulse masses)
  propsd[2]  = Rf_asReal(pv_sd_pulse_mass);  // Proposal variance for SD of RE masses (formerly pv_pulse_mass)
  propsd[3]  = Rf_asReal(pv_sd_pulse_width); // Proposal variance for SD of RE masses
  propsd[6]  = Rf_asReal(pv_pulse_location); 
  propsd[0]  = Rf_asReal(pv_baseline); 
  propsd[1]  = Rf_asReal(pv_halflife);
  propsd[7]  = Rf_asReal(pv_mean_pulse_mass);  // sd of the proposed fixed effect mass distr
  propsd[8]  = Rf_asReal(pv_mean_pulse_width); // sd of the proposed fixed effect width distr
  propsd[9]  = Rf_asReal(pv_etamass);  // sd of the proposed eta mass distribution
  propsd[10] = Rf_asReal(pv_etawidth); // sd of the proposed eta width distribution

  //Rprintf("proposal variance for location is %f\n", propsd[6]);
  //Rprintf("\nProposal sd indiv mass: %f", propsd[4]);
  //Rprintf("\nProposal sd indiv width: %f", propsd[5]);
  //Rprintf("\nProposal sd sd mass: %f", propsd[2]);
  //Rprintf("\nProposal sd sd width: %f", propsd[3]);
  //Rprintf("\nProposal sd location: %f", propsd[6]);
  //Rprintf("\nProposal sd baseline: %f", propsd[0]);
  //Rprintf("\nProposal sd halflife: %f", propsd[1]);
  //Rprintf("\nProposal sd mean mass: %f", propsd[7]);
  //Rprintf("\nProposal sd mean width: %f", propsd[8]);
  //Rprintf("\nProposal sd eta mass: %f", propsd[9]);
  //Rprintf("\nProposal sd eta width: %f", propsd[10]);


  // Set up pulse linklist ------------------
  Node_type *list;      
  list = initialize_node();


  // R lists for chains ------------------
  // Common parameters -- 1 record per saved iteration
  SEXP common1 = Rf_protect(Rf_allocMatrix(REALSXP, out_len, 9));
  // Pulse-specific parameters -- 1 record per pulse per iteration
  //   list object defined here, each iteration will differ in length, so the
  //   entries are defined in linklistv3
  SEXP pulse_chains = Rf_protect(Rf_allocVector(VECSXP, out_len));


  mcmc(list, parms, ts, iters, nobs, nthin, nburnin, strauss, verbose, priors, common1,
       pulse_chains, propsd);


  // Column names for common1 chain --------
  SEXP dim, dimnames;
  // create 'dim' attribute (dimension of matrix)
  Rf_protect(dim = Rf_allocVector(INTSXP, 2));
  INTEGER(dim)[0] = out_len; 
  INTEGER(dim)[1] = 9;
  Rf_setAttrib(common1, R_DimSymbol, dim);
  // create column names for common parms matrix
  SEXP common_names = Rf_protect(Rf_allocVector(STRSXP, 9));
  SET_STRING_ELT(common_names, 0, Rf_mkChar("iteration"));
  SET_STRING_ELT(common_names, 1, Rf_mkChar("num_pulses"));
  SET_STRING_ELT(common_names, 2, Rf_mkChar("baseline"));
  SET_STRING_ELT(common_names, 3, Rf_mkChar("mean_pulse_mass"));
  SET_STRING_ELT(common_names, 4, Rf_mkChar("mean_pulse_width"));
  SET_STRING_ELT(common_names, 5, Rf_mkChar("halflife"));
  SET_STRING_ELT(common_names, 6, Rf_mkChar("model_error"));
  SET_STRING_ELT(common_names, 7, Rf_mkChar("sd_mass"));
  SET_STRING_ELT(common_names, 8, Rf_mkChar("sd_widths"));
  // assign names to vector
  Rf_protect(dimnames = Rf_allocVector(VECSXP, 2));
  // assign names vector to columns (1) names attribute, leaving rows (0) NULL;
  SET_VECTOR_ELT(dimnames, 1, common_names);
  // assign names vector to common1 chain
  Rf_setAttrib(common1, R_DimNamesSymbol, dimnames);


  // Combine chains into list for output -------------
  SEXP chains = Rf_protect(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(chains, 0, common1);
  SET_VECTOR_ELT(chains, 1, pulse_chains);


  // Free memory ---------------------------
  deallocate_data(ts, nobs);
  destroy_list(list);
  free(priors->fe_mean);
  free(priors->fe_variance);
  free(priors->re_sdmax);
  free(priors);
  //free(parms->re_sd);
  free(parms);
  Rf_unprotect(6);

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
    data[i][0] = REAL(getListElement(indata, "time"))[i];
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

