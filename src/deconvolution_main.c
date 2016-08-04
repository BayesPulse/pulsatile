//-----------------------------------------------------------------------------
//
// MAIN BDMCMC PROGRAM
// 
// FILE: deconvolution_main.c
// AUTHOR: Matt Mulvahill 
//         (and many before me: T Johnson, N Carlson, K Horton, Karen )
//
// DESCRIPTION: 
//   Contains MAIN function for bdmcmc deconvolution analysis.  
//   
//   This version of the program can use one of two types of BD algorithms:
// 
//     1) Standard BDMCMC as described by Stephens 2000, using an 
//        order-statistic prior for pulse location and poisson prior for pulse
//        count. To use this version, set parms->gamma = -1 in the argument
//        file (parms->gamma is read in as priorgamma). The order-statistic to
//        use is set by 'mmm'
//     2) Spatial BDMCMC as described in Moller's 'Statistical Inference and
//        Simulation for Spatial Point Processes.' In our algorithm, the prior
//        for pulse count and location is set as a repulsive Strauss process. 
//        To use the Strauss prior, set parms->gamma >= 0 and <= 1 and
//        parms->range to the desired number of minutes.  A gamma value of 0
//        is a Hard Core process, where no other pulse can be within +/-range 
//        of another pulse.
//
// Notes:
//   FE and RE stand for Fixed and Random Effects
// 
//   Note that a gamma > 1 will throw an error.
// 
// SUBROUTINES: 
//   main()
//   
// GLOBAL VARIABLE DEFINITIONS:
//   fitstart - The first time in hours that a pulse may occur
//   fitend   - The last time in hours that a pulse may occur
//   mmm      - Order statistic used for distribution of pulse locations.
//              This is assigned in deconvolution_main.c and is typically 3.
// 
//-----------------------------------------------------------------------------

#define R_NO_REMAP 
#include <R.h>
#include <Rinternals.h>
#include "deconvolution_main.h"
#include "format_data.h"
#include "hash.h"
#include "linklistv3.h"
#include "birthdeath_strauss.h"
#define EPS 1.0e-12 // error for testing equality with 0

// Global variables --------------------
double M; // Precision used in RNGs. Set to 32 or 64 based on operating system.
int mmm = 3;     // Order statistic to use for order stat prior
double fitstart; // First time a pulse can occur (10 min. increments).
double fitend;   // Last time a pulse can occur (10 min. increments).




//
// Primary analysis function  ('main')
// 
// Removed SEXP thin and SEXP seeds from 
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


  // Declarations ---------------------
  int i;                          // Generic counter
  int *N;                         // Number of obs in datafile, set by read_data_file()
  int iter = INTEGER(iterations); // number of iterations to run mcmc (args file)
  int a;                          // Generic counter
  int placeholder = 0;            // For return value on fscanf's (ignoring them)
  //unsigned long *seed;          // Pointer to seeds (args file)
  double **ts;                    // Pointer to Nx2 matrix containing data
  double propvar[7];              // Array of proposal variances for MH algorithms
  double time[9];                 // Starting values for pulse locations
  double mass[9];                 // Starting values for pulse masses
  double width[9];                // Starting values for pulse widths

  Common_parms *parms;            // Common parameters data structure
  Priors *priors;                 // Prior parameters data structure
  Node_type *list;                // Linked list of nodes/pulses
  Node_type *new_node;            // New node/pulse

  // Rprintf(parm, "%d %d %d %lf %lf %lf\n", 
  //        i/NN, num_node2, num_node, new_node->time, 
  //        new_node->theta[0], new_node->theta[1]);
  // num_node++;
  // new_node = new_node->succ;
  //
  // Rprintf(common, "%d %lf %lf %lf %lf %lf %lf %lf \n", 
  //         num_node2, parms->md[0], parms->theta[0], parms->theta[1],
  //         parms->md[1], parms->sigma, parms->re_sd[0], parms->re_sd[1]);
  
  // R lists for returning chains
  // Common parameters -- 1 record per saved iteration
  SEXP common1 = Rf_protect(Rf_allocMatrix(REALSXP, iter/NN, 8));
  // Pulse-specific parameters -- 1 record per pulse per iteration
  //   list object defined here, each iteration will differ in lenght, so the
  //   entries are defined in linklistv3
  SEXP parm1   = PROTECT(Rf_allocVector(VECSXP, iter/NN));
  //char common1[100];   // File path to save common parameters
  //char parm1[100];     // File path to save pulse-specific parameters

  // Start main function
  for (a = 0; a < 1; a++) {

    //--------------------------------------------
    // Read arguments file and print to log
    //--------------------------------------------
    
    // Allocate memory for seed array
    //seed  = (unsigned long *)calloc(3, sizeof(unsigned long));


    // -------------------------------------------
    // Read in the hormonal time series 
    // -------------------------------------------
    ts = convert_data(indata);


    //--------------------------------------------
    // Allocate memory and initialize values 
    //--------------------------------------------
    
    // Precision used in KISS RNG --------------
    //M = exp(-64 * log(2.0));

    // Range of real time for model fit ---------
    fitend   = ts[*N - 1][0] + ts[0][0] * 2; // search 2 units farther in time
    fitstart = -ts[0][0] * 4;                // search 4 units in the past

    // Set up priors structure ------------------
    priors                 = (Priors *)calloc(1, sizeof(Priors));
    priors->re_sdmax       = (double *)calloc(2, sizeof(double));
    priors->fe_variance    = (double *)calloc(2, sizeof(double));
    priors->fe_mean        = (double *)calloc(2, sizeof(double));
    priors->meanbh[0]      = REAL(prior_baseline_mean);    // priormub;
    priors->meanbh[1]      = REAL(prior_halflife_mean);    // priormuh;
    priors->varbh[0]       = REAL(prior_baseline_var);     // priorvarb;
    priors->varbh[1]       = REAL(prior_halflife_var);     // priorvarh;
    priors->fe_mean[0]     = REAL(prior_pulse_mass_mean);  // priormu1;
    priors->fe_mean[1]     = REAL(prior_pulse_width_mean); // priormu2;
    priors->fe_variance[0] = REAL(prior_pulse_mass_var);   // priorvar1;
    priors->fe_variance[1] = REAL(prior_pulse_width_var);  // priorvar2;
    priors->re_sdmax[0]    = REAL(prior_max_sd_mass);      // priora1;
    priors->re_sdmax[1]    = REAL(prior_max_sd_width);     // priora2;
    priors->err_alpha      = REAL(prior_error_alpha);      // prioralpha;
    priors->err_beta       = REAL(prior_error_beta);       // priorbeta;
    priors->gamma          = REAL(prior_pulse_location_gamma);   // priorgamma;
    priors->range          = REAL(prior_pulse_location_range);   // priorrange;

    // Set up parms structure -------------------
    parms           = (Common_parms *)calloc(1, sizeof(Common_parms));
    parms->re_sd    = (double *)calloc(2, sizeof(double));
    parms->nprior   = REAL(prior_pulse_location_count); // priorr;
    parms->theta[0] = REAL(sv_pulse_mass_mean);   // svmu1;
    parms->theta[1] = REAL(sv_pulse_width_mean);  // svmu2;
    parms->md[0]    = REAL(sv_baseline_mean);     // svbase;
    parms->md[1]    = REAL(sv_halflife_mean);     // svhalf;
    parms->sigma    = REAL(sv_error_var);         // svevar; // error variance
    parms->lsigma   = REAL(log(parms->sigma));    // log of error variance
    parms->re_sd[0] = REAL(sv_pulse_mass_sd);  // svsig1; 
    parms->re_sd[1] = REAL(sv_pulse_width_sd); // svsig2; 
    // note: other re_sd parms (max and pv) are entered and used on the natural
    // scale, so for consistency, we enter these on the same scale.  Functions
    // other than draw_re_sd, use re_sd[j] on the log scale so we log xform them
    // here.                                                                   
    parms->decay    = log(2) / parms->md[1]; // calculate decay rate from halflife


    //
    // Initialize pulse values 
    //   NOTE: Possibly add option to pulse_spec for setting starting values/num
    //   pulses
    //
    list = initialize_node();
    mass[0] = 12.32657558;
    mass[1] = 8.335093817;
    mass[2] = 6.169354435;
    mass[3] = 17.80691618;
    mass[4] = 17.81203064;
    mass[5] = 4.180821883;
    mass[6] = 9.385968508;
    mass[7] = 3.539176735;
    mass[8] = 8.860152668;
    width[0] = 26.02652032;
    width[1] = 37.383667;
    width[2] = 37.90962662;
    width[3] = 31.45537546;
    width[4] = 52.64992824;
    width[5] = 42.00314896;
    width[6] = 11.7825812;
    width[7] = 21.11221475;
    width[8] = 100.6173752;
    time[0] = 139.2555248;
    time[1] = 293.9643442;
    time[2] = 400.3647507;
    time[3] = 654.7396597;
    time[4] = 788.742805;
    time[5] = 918.9406528;
    time[6] = 1156.007127;
    time[7] = 1335.413248;
    time[8] = 1408.372033;

    // Initialize nodes and insert into linkedlist
    for (i = 0; i < 9; i++) {

      new_node           = initialize_node();
      new_node->time     = time[i];
      new_node->theta[0] = mass[i];
      new_node->theta[1] = width[i];
      new_node->mean_contrib = (double *)calloc(*N, sizeof(double));
      mean_contribution(new_node, ts, parms, *N);
      insert_node(new_node, list);

    }



    //
    // Run MCMC
    //
    GetRNGstate();
    mcmc(list, parms, ts, iter, *N, priors, common1, parm1, propvar);
    PutRNGstate();


    //
    // Deallocate resources 
    //
    destroy_list(list);
    free(seed);
    for (i = 0; i <* N; i++) {
      free(ts[i]);
    }
    free(ts);
    free(N);
    free(priors->fe_mean);
    free(priors->fe_variance);
    free(parms->re_sd);
    free(priors->re_sdmax);
    free(priors);
    free(parms);

  }

  // Exit with success code
  exit(EXIT_SUCCESS);

}



//
// convert data from SEXP to C array
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


//-----------------------------------------------------------------------------
// End of file
//-----------------------------------------------------------------------------

