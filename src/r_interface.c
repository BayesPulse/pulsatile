//
// Primary interface between R and C code
//


//#define R_NO_REMAP 
//#include <R.h>
//#include <Rinternals.h>
// For boolean variables, consider: #include <stdbool.h>
#include <stdlib.h>
#include <stdbool.h>
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

  fitend   = ts[nobs - 1][0] + ((ts[1][0] - ts[0][0]) * 2); // search 2 units farther in time
  fitstart = -(ts[1][0] - ts[0][0]) * 4;                    // search 4 units in the past

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
  priors->err_beta       = 1./Rf_asReal(prior_error_beta);    // priorbeta -- input as rate, convert to scale for Rf_rgamma();
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

  // NOTE: FOR TESTING ONLY
  bool TEST = true;

  // If testing, use hard-coded pulse data
  if (TEST) {
    list = test_pulses(list, ts, parms, nobs);
  } else {
    list = initialize_node();
  }


  // R lists for chains ------------------
  // Common parameters -- 1 record per saved iteration
  SEXP common1 = Rf_protect(Rf_allocMatrix(REALSXP, out_len, 9));
  // Pulse-specific parameters -- 1 record per pulse per iteration
  //   list object defined here, each iteration will differ in length, so the
  //   entries are defined in linklistv3
  SEXP pulse_chains = Rf_protect(Rf_allocVector(VECSXP, out_len));

  double tmp = likelihood(list, ts, parms, nobs, list);
  Rprintf("Initial likelihood is: %f\n", tmp);


  mcmc(list, parms, ts, iters, nobs, nthin, nburnin, strauss, verbose, priors,
       common1, pulse_chains, propsd);


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

//Node_type * test_data(Node_type * list, int nobs) {
//
//  Rcpp::NumericVector thistime(144);
//  for (int i = 0; i < thistime.size(); i++)  thistime(i) = (i + 1) * 10;
//      Rcpp::NumericVector conc =
//      { 2.681436, 3.619304, 4.583596, 5.391640, 4.602580, 5.837577, 5.093759,
//        4.460362, 4.228568, 3.749578, 3.756091, 3.794117, 3.458053, 3.614573,
//        2.979374, 3.235577, 3.379498, 3.083957, 3.299427, 3.467285, 2.969050,
//        3.401228, 3.146146, 2.934315, 3.041346, 3.161088, 2.740992, 2.869972,
//        3.236580, 5.192121, 5.752548, 5.496779, 5.526355, 4.945300, 4.936505,
//        7.312601, 7.128395, 6.801745, 7.635689, 6.246128, 6.032690, 4.702559,
//        5.039149, 4.219767, 4.527492, 4.090115, 3.455401, 3.187411, 4.271770,
//        6.863361, 5.953694, 6.456755, 6.778942, 4.981977, 4.830248, 4.336557,
//        4.296211, 4.237477, 3.427465, 3.665419, 2.978837, 3.409569, 2.641762,
//        3.689739, 3.180302, 3.497940, 3.097911, 3.205925, 4.962560, 5.559928,
//        6.377035, 4.974151, 5.094906, 4.864500, 4.852309, 5.315221, 5.491787,
//        5.379786, 4.845070, 4.835242, 4.523089, 4.211985, 4.741966, 6.410714,
//        7.966407, 7.233676, 6.293541, 5.807553, 5.626408, 4.685379, 4.976104,
//        4.923761, 5.616314, 4.954546, 4.316296, 4.449381, 4.035612, 5.933037,
//        5.464214, 5.145751, 5.200191, 4.553076, 4.429967, 3.915830, 3.962575,
//        3.418965, 3.334863, 3.174500, 3.409328, 2.822615, 3.298277, 2.421233,
//        3.413683, 2.850547, 3.115562, 2.713616, 2.941980, 2.887866, 2.980766,
//        3.627824, 4.625605, 4.468451, 4.815758, 3.985436, 3.463471, 3.682286,
//        3.536958, 3.563942, 3.552810, 3.751709, 2.933170, 4.234158, 4.716445,
//        4.043727, 4.320064, 3.972299, 3.851225, 4.050221, 3.195143, 3.168399,
//        3.011654, 2.721384, 3.279211, 3.079000 };
//
//  // Declare variables
//  int i = 0;
//  double **data; // Nx2 matrix of time and log-concentration
//
//  // Allocate memory for matrix of data
//  data = (double **)calloc(nrow, sizeof(double *));
//  for (i = 0; i < nrow; i++) {
//    data[i] = (double *)calloc(2, sizeof(double));
//  }
//
//  // Read the time series into memory                          
//  for (i = 0; i < nrow; i++) {
//    data[i][0] = REAL(getListElement(indata, "time"))[i];
//    data[i][1] = REAL(getListElement(indata, "concentration"))[i];
//    data[i][1] = log(data[i][1]);
//  }
//
//	return data;
//}


Node_type * test_pulses(Node_type * list, double **ts, Common_parms *parms, 
                        int nobs) {

  list = initialize_node();

  // Initialize pulse values -------------------
  double time[12];      // Starting values for pulse locations
  double mass[12];      // Starting values for pulse masses
  double width[12];     // Starting values for pulse widths
  double eta_mass[12];  // Starting values for eta masses
  double eta_width[12]; // Starting values for eta widths
  mass[0]       = 4.2407384;
  mass[1]       = 0.7655571;
  mass[2]       = 4.0825812;
  mass[3]       = 4.7655685;
  mass[4]       = 4.7903024;
  mass[5]       = 4.0835315;
  mass[6]       = 2.3587096;
  mass[7]       = 5.3756180;
  mass[8]       = 1.2616455;
  mass[9]       = 2.7332376;
  mass[10]      = 2.4361395;
  mass[11]      = 2.1984437;
  width[0]      = 88.111121;
  width[1]      = 110.366726;
  width[2]      = 36.036933;
  width[3]      = 49.107879;
  width[4]      = 20.283757;
  width[5]      = 36.461844;
  width[6]      = 57.384845;
  width[7]      = 57.401967;
  width[8]      = 1.929962;
  width[9]      = 15.429365;
  width[10]     = 139.865751;
  width[11]     = 43.961197;
  time[0]       = 26.54152;
  time[1]       = 174.63993;
  time[2]       = 298.62117;
  time[3]       = 360.55329;
  time[4]       = 494.61155;
  time[5]       = 689.09242;
  time[6]       = 763.89017;
  time[7]       = 839.80027;
  time[8]       = 925.80251;
  time[9]       = 975.47320;
  time[10]      = 1199.00866;
  time[11]      = 1322.82471;
  eta_mass[0]   = 0.8486756;
  eta_mass[1]   = .7108269;
  eta_mass[2]   = .9045537;
  eta_mass[3]   = .4203106;
  eta_mass[4]   = .8413771;
  eta_mass[5]   = .4562658;
  eta_mass[6]   = .2856221;
  eta_mass[7]   = .2336739;
  eta_mass[8]   = .6728940;
  eta_mass[9]   = .9409146;
  eta_mass[10]  = .6038081;
  eta_mass[11]  = .6312320;
  eta_width[0]  = 0.6803521;
  eta_width[1]  = 0.6169701;
  eta_width[2]  = 0.7486873;
  eta_width[3]  = 0.4655990;
  eta_width[4]  = 1.3311321;
  eta_width[5]  = 0.5488085;
  eta_width[6]  = 0.6290294;
  eta_width[7]  = 1.4267846;
  eta_width[8]  = 0.6693473;
  eta_width[9]  = 0.5769079;
  eta_width[10] = 0.5871716;
  eta_width[11] = 0.8507986;



  // Initialize nodes and insert into linkedlist
  for (int i = 0; i < 12; i++) {

    Node_type *new_node;
    new_node               = initialize_node();
    new_node->time         = time[i];
    new_node->theta[0]     = mass[i];
    new_node->theta[1]     = width[i];
    new_node->eta[0]       = eta_mass[i];
    new_node->eta[1]       = eta_width[i];
    new_node->mean_contrib = (double *)calloc(nobs, sizeof(double));
    mean_contribution(new_node, ts, parms, nobs);
    insert_node(new_node, list);

  }

  return list;

}


//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

