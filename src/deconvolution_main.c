///-----------------------------------------------------------------------------
///
/// MAIN BDMCMC PROGRAM
/// 
/// FILE: deconvolution_main.c
/// AUTHOR: Matt Mulvahill 
///         (and many before me: T Johnson, N Carlson, K Horton, Karen )
///
/// DESCRIPTION: 
///   Contains MAIN function for bdmcmc deconvolution analysis.  
///   
///   This version of the program can use one of two types of BD algorithms:
/// 
///     1) Standard BDMCMC as described by Stephens 2000, using an 
///        order-statistic prior for pulse location and poisson prior for pulse
///        count. To use this version, set parms->gamma = -1 in the argument
///        file (parms->gamma is read in as priorgamma). The order-statistic to
///        use is set by 'mmm'
///     2) Spatial BDMCMC as described in Moller's 'Statistical Inference and
///        Simulation for Spatial Point Processes.' In our algorithm, the prior
///        for pulse count and location is set as a repulsive Strauss process. 
///        To use the Strauss prior, set parms->gamma >= 0 and <= 1 and
///        parms->range to the desired number of minutes.  A gamma value of 0
///        is a Hard Core process, where no other pulse can be within +/-range 
///        of another pulse.
///
/// Notes:
///   FE and RE stand for Fixed and Random Effects
/// 
///   Note that a gamma > 1 will throw an error.
/// 
/// SUBROUTINES: 
///   main()
///   
/// GLOBAL VARIABLE DEFINITIONS:
///   fitstart - The first time in hours that a pulse may occur
///   fitend   - The last time in hours that a pulse may occur
///   mmm      - Order statistic used for distribution of pulse locations.
///              This is assigned in deconvolution_main.c and is typically 3.
/// 
///-----------------------------------------------------------------------------

#include <R.h>
#include <Rinternals.h>
#include "deconvolution_main.h"
#include "format_data.h"
#include "hash.h"
#include "linklistv3.h"
#include "birthdeath_strauss.h"
//#include <gperftools/profiler.h> // for profiling code (time spent in each function)

#define EPS 1.0e-12 // error for testing equality with 0

// Global variables --------------------
double M;        // Precision used in RNGs. Set to 32 or 64 based on operating system.
int mmm = 3;     // Order statistic to use, set regardless of using orderstat or strauss
double fitstart; // First time a pulse can occur (10 min. increments).
double fitend;   // Last time a pulse can occur (10 min. increments).


///-----------------------------------------------------------------------------
/// 
///  Call function for main program decon
/// 
/// int main(int argc, char *argv[]) {
///-----------------------------------------------------------------------------
//decon_call


///-----------------------------------------------------------------------------
/// 
///  Main program starts here
/// 
/// int main(int argc, char *argv[]) {
///-----------------------------------------------------------------------------
SEXP decon(SEXP dataset, 
           SEXP thin, 
           SEXP prior_parms, 
           SEXP starting_vals, 
           SEXP proposal_vars, 
           SEXP seeds) {


  // Declarations ---------------------
  int i;               // Generic counter
  int *N;              // Number of obs in datafile, set by read_data_file()
  int iter;            // number of iterations to run mcmc (args file)
  int a;               // Generic counter
  int placeholder = 0; // For return value on fscanf's (ignoring them)
  //char datafile[90];   // File path to datafile
  char common1[100];   // File path to save common parameters
  char parm1[100];     // File path to save pulse-specific parameters
  unsigned long *seed; // Pointer to seeds (args file)
  double **ts;         // Pointer to Nx2 matrix containing data
  double propvar[7];   // Array of proposal variances for MH algorithms
  double priormu1;     // Prior mean pulse mass FE
  double priorvar1;    // Prior variance of mean pulse mass FE
  double priormu2;     // Prior mean pulse width FE
  double priorvar2;    // Prior variance of mean pulse width FE
  double priormub;     // Prior mean baseline concentration
  double priorvarb;    // Prior variance of mean baseline concentration
  double priormuh;     // Prior mean half-life
  double priorvarh;    // Prior variance of mean half-life
  double prioralpha;   // Prior alpha on inverse gamma for error
  double priorbeta;    // Prior beta on inverse gamma for error
  double priora1;      // Maximum SD for Uniform prior on variance of mass RE
  double priora2;      // Maximum SD for Uniform prior on variance of width RE
  double priorr;       // Prior num. pulses/pulse intensity on Poisson/Strauss
  double priorgamma;   // Prior for gamma in Strauss
  double priorrange;   // Prior for range in Strauss
  double svmu1;        // Starting value for mean pulse mass FE
  double svmu2;        // Starting value for mean pulse width FE
  double svbase;       // Starting value for baseline conc
  double svhalf;       // Starting value for half-life
  double svevar;       // Starting value for error variance
  double svsig1;       // Starting value for RE mass standard deviation
  double svsig2;       // Starting value for RE width standard deviation
  double time[9];      // Starting values for pulse locations
  double mass[9];      // Starting values for pulse masses
  double width[9];     // Starting values for pulse widths

  //FILE *finput;        // Arguments file
  Common_parms *parms; // Common parameters data structure
  Priors *priors;      // Prior parameters data structure
  Node_type *list;     // Linked list of nodes/pulses
  Node_type *new_node; // New node/pulse


  // Start main function
  for (a=0;a<1;a++) {

    // Print decon version based on git commits 
    // remove these lines if you aren't  using git to manage the source code
    //Rprintf("Welcome to decon.\n");
    //Rprintf("You are using version %s\n", VERSION);
    //Rprintf("\n\n");


    //--------------------------------------------
    // Read arguments file and print to log
    //--------------------------------------------
    
    // Allocate memory for seed array
    seed  = (unsigned long *)calloc(3, sizeof(unsigned long));

    // Open arguments file for reading
    //finput = fopen(argv[1], "r");
    //if (finput == NULL) {
    //  perror("Argument error");
    //  Rprintf(stderr, "Value of argument: %s\n", argv[1]);
    //  exit(EXIT_FAILURE);
    //}


    // Read arguments file line by line
    //placeholder = fscanf(finput, "%s \n", datafile);
    //placeholder = fscanf(finput, "%s \n", common1);
    //placeholder = fscanf(finput, "%s \n", parm1);
    //placeholder = fscanf(finput, "%lu %lu %lu\n", &seed[0], &seed[1], &seed[2]);

    //placeholder = fscanf(finput, "%d \n", &iter);
    //placeholder = fscanf(finput, "%lf %lf\n", &priormu1, &priorvar1);
    //placeholder = fscanf(finput, "%lf %lf\n", &priormu2, &priorvar2);
    //placeholder = fscanf(finput, "%lf %lf\n", &priormub, &priorvarb);
    //placeholder = fscanf(finput, "%lf %lf\n", &priormuh, &priorvarh);
    //placeholder = fscanf(finput, "%lf %lf\n", &prioralpha, &priorbeta);
    //placeholder = fscanf(finput, "%lf %lf\n", &priorgamma, &priorrange);
    //placeholder = fscanf(finput, "%lf %lf\n", &priora1, &priora2);
    //placeholder = fscanf(finput, "%lf\n",  &priorr);
    //placeholder = fscanf(finput, "%lf %lf\n", &svmu1, &svmu2);
    //placeholder = fscanf(finput, "%lf %lf\n", &svbase, &svhalf);
    //placeholder = fscanf(finput, "%lf\n", &svevar);
    //placeholder = fscanf(finput, "%lf %lf\n", &svsig1, &svsig2);
    //placeholder = fscanf(finput, "%lf %lf\n", &propvar[0], &propvar[1]);
    //placeholder = fscanf(finput, "%lf %lf\n", &propvar[2], &propvar[3]);
    //placeholder = fscanf(finput, "%lf %lf %lf\n", &propvar[4], &propvar[5], &propvar[6]);
    //fclose(finput);

    // Print arguments to log for confirming what was used
    //Rprintf("Argument file used:\n");
    //Rprintf("%s \n", argv[1]);
    //Rprintf("Contents of argument file:\n");
    //Rprintf("%s \n", datafile);
    //Rprintf("%s \n", common1);
    //Rprintf("%s \n", parm1);
    //Rprintf("%lu %lu %lu\n", seed[0], seed[1], seed[2]);
    //Rprintf("%d \n", iter);
    //Rprintf("%lf %lf\n",  priormu1,  priorvar1);
    //Rprintf("%lf %lf\n",  priormu2,  priorvar2);
    //Rprintf("%lf %lf\n",  priormub,  priorvarb);
    //Rprintf("%lf %lf\n",  priormuh,  priorvarh);
    //Rprintf("%lf %lf\n",  prioralpha,  priorbeta);
    //Rprintf("%lf %lf\n",  priorgamma, priorrange);
    //Rprintf("%lf %lf\n",  priora1,  priora2);
    //Rprintf("%lf\n",  priorr);
    //Rprintf("%lf %lf\n", svmu1, svmu2);
    //Rprintf("%lf %lf\n", svbase, svhalf);
    //Rprintf("%lf\n", svevar);
    //Rprintf("%lf %lf\n", svsig1, svsig2);
    //Rprintf("%lf %lf\n", propvar[0], propvar[1]);
    //Rprintf("%lf %lf\n", propvar[2], propvar[3]);
    //Rprintf("%lf %lf %lf\n\n\n", propvar[4], propvar[5], propvar[6]);



    // -------------------------------------------
    // Read in the hormonal time series 
    // -------------------------------------------
    //N   = (int *)calloc(1, sizeof(int));
    //ts  = read_data_file(datafile, N);


    // -------------------------------------------
    // Print which prior will be used
    // -------------------------------------------
    if (priorgamma < -0.001) {
      Rprintf("Using the %d-order statistic prior.\n", mmm);

    } else if (priorgamma > 1.001) {

      Rprintf("Invalid value used for priorgamma in Strauss.\n");
      Rprintf("priorgamma = %f", priorgamma); 
      Rprintf("Error occurred prior to running mcmc\n");
      exit(EXIT_FAILURE);
      
    } else {

      Rprintf("Using the Strauss prior with gamma = %f\n", priorgamma);

      if (fabs(priorgamma) < EPS && priorgamma != 0) {
        Rprintf("WARNING: Hard Core attempted, but priorgamma NOT");
        Rprintf("exactly equal to 0!!\n");

      }

    }


    //--------------------------------------------
    // Allocate memory and initialize values 
    //--------------------------------------------
    
    // Precision used in KISS RNG --------------
    M = exp(-64 * log(2.0));

    // Range of real time for model fit ---------
    fitend   = ts[*N - 1][0] + ts[0][0] * 2; // search 2 units farther in time
    fitstart = -ts[0][0] * 4;                // search 4 units in the past

    // Set up priors structure ------------------
    priors                 = (Priors *)calloc(1, sizeof(Priors));
    priors->re_sdmax       = (double *)calloc(2, sizeof(double));
    priors->fe_variance    = (double *)calloc(2, sizeof(double));
    priors->fe_mean        = (double *)calloc(2, sizeof(double));
    priors->meanbh[0]      = priormub;
    priors->meanbh[1]      = priormuh;
    priors->varbh[0]       = priorvarb;
    priors->varbh[1]       = priorvarh;
    priors->fe_mean[0]     = priormu1;
    priors->fe_mean[1]     = priormu2;
    priors->fe_variance[0] = priorvar1;
    priors->fe_variance[1] = priorvar2;
    priors->re_sdmax[0]    = priora1;
    priors->re_sdmax[1]    = priora2;
    priors->err_alpha      = prioralpha;
    priors->err_beta       = priorbeta;
    priors->gamma          = priorgamma;
    priors->range          = priorrange;

    // Set up parms structure -------------------
    parms           = (Common_parms *)calloc(1, sizeof(Common_parms));
    parms->re_sd    = (double *)calloc(2, sizeof(double));
    parms->nprior   = priorr;
    parms->theta[0] = svmu1;
    parms->theta[1] = svmu2;
    parms->md[0]    = svbase;
    parms->md[1]    = svhalf;
    parms->sigma    = svevar;  // error variance
    parms->lsigma   = log(parms->sigma); // log of error variance
    parms->re_sd[0] = svsig1; // note: other re_sd parms (max and pv) are 
    parms->re_sd[1] = svsig2; // entered and used on the natural scale, so 
                                   // for consistency, we enter these on the
                                   // same scale.  Functions other than
                                   // draw_re_sd, use re_sd[j] on the log scale
                                   // so we log xform them here.
    parms->decay    = log(2) / parms->md[1]; // calculate decay rate from halflife

    // Initialize pulse values -------------------
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



    //-------------------------------------------
    // Run MCMC
    //-------------------------------------------
    //ProfilerStart("decon.prof");
    mcmc(list, parms, ts, iter, *N, priors, seed, common1, parm1, propvar);
    //ProfilerStop();


    //-------------------------------------------
    // Deallocate resources 
    //-------------------------------------------
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





//-----------------------------------------------------------------------------
// End of file
//-----------------------------------------------------------------------------

