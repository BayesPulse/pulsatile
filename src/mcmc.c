//-----------------------------------------------------------------------------
//
// FILE: mcmc.c
// MAINTAINER: Matt Mulvahill 
//         (and many before me: T Johnson, N Carlson, K Horton, Karen )
//
// DESCRIPTION: 
//   Contains birth-death algorithm, expanded to spatial birth-death with the
//   option to use the order-statistic or strauss prior on pulse location.  To
//   use the order statistic, set strauss=0; strauss=1 for the strauss prior. 
// 
//------------------------------------------------------------------------------

// Include needed header files
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "mcmc.h"
#include "r_interface.h"
#include "pulse_node.h"
#include "birth_death.h"
#include "calculations.h"

// Minumum precision, used instead of equalities in some places
#define EPS 1.0e-42

// Global variables
extern double fitstart; // First time in 10min increments a pulse can occur
extern double fitend;   // Last time in 10min increments a pulse can occur




//-----------------------------------------------------------------------------
//
// mcmc()
//   This function runs the BDMCMC process
//
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that exist
//     Common_parms *parms - the current values of the common parameters 
//     double **ts         - this is the matrix of observed data (a column
//                           of times and a column of log(concentration) 
//     long iter           - the number of iterations to run
//     int N               - the number of observations in **ts
//     int nthin           - Keep every nthin sample
//     int burnin          - Number of iterations to remove for burnin
//     int strauss         - indicator for whether to use the strauss (=1) or
//                           order-stat (=0) pulse location prior 
//     Priors *priors      - the parameters of the prior distributions
//     char *file1         - the output file name for common parameters
//     char *file2         - the output file name for pulse specific
//                           parameters
//     double **pmd_var    - proposal variance-covariance matrix for
//                           baseline and halflife
//
//   RETURNS: 
//     None                - all updates are made internally
// 
//-----------------------------------------------------------------------------
void mcmc(Node_type *list, 
          Common_parms *parms, 
          double **ts, 
          long iter,
          int N,
          int nthin,
          int burnin,
          int strauss,
          int verbose,
          Priors *priors,
          SEXP common, //char *file1, 
          SEXP pulse_chains, //char *file2, 
          double propsd[]) {

  //
  // Declarations 
  //
  GetRNGstate();
  int i;               // Generic counter
  int j;               // Generic counter
  int k;               // Generic counter
  int l;               // Generic counter
  int num_node;        // Counter of num_node for loops
  int num_node2;       // Number of pulses
  int out_len = (iter-burnin) / nthin;
  int out_i;
  int NNN = 5000;      // Output ever NNth iteration to STDOUT
  double sdrem;        // Proposal SD for RE masses
  double sdrew;        // Proposal SD for RE widths
  double sdmv;         // Proposal SD for SD of RE masses
  double sdwv;         // Proposal SD for SD of RE widths
  double sdfem;        // Proposal SD for FE masses
  double sdfew;        // Proposal SD for FE widths
  double sdt;          // Proposal variance for individual pulse locations
  double ssq;          // Sum of squared differences between log(concentration) and expected value
  double sdetam;
  double sdetaw;
  double *likeli;      // Value of likelihood
  double **pmd_var;    // Var-cov matrix for baseline and half-life
  double **pmd_vch;    // Cholesky decomposed var-cov matrix for baseline/half-life
  Node_type *new_node; // Node structure for new pulse

  //
  // Counters for acceptences (a) and draws (n)
  //
  long adelta = 0; // halflife MH
  long ndelta = 0; // halflife MH
  long atime  = 0; // pulse location MH
  long ntime  = 0; // pulse location MH
  long arem   = 0; // RE pulse mass
  long nrem   = 0; // RE pulse mass
  long arew   = 0; // RE pulse width
  long nrew   = 0; // RE pulse width
  long arevm  = 0; // RE pulse mass SD
  long nrevm  = 0; // RE pulse mass SD
  long arevw  = 0; // RE pulse width SD
  long nrevw  = 0; // RE pulse width SD
  long afem   = 0; // Fixed Effect pulse mass
  long nfem   = 0; // Fixed Effect pulse mass
  long afew   = 0; // Fixed Effect pulse width
  long nfew   = 0; // Fixed Effect pulse width
  long aetam  = 0; // pulse mass eta 
  long aetaw  = 0; // pulse width eta
  long netam  = 0; // pulse mass eta 
  long netaw  = 0; // pulse width eta
  long *adelta_ptr = &adelta; // halflife MH
  long *ndelta_ptr = &ndelta; // halflife MH
  long *atime_ptr  = &atime ; // pulse location MH
  long *ntime_ptr  = &ntime ; // pulse location MH
  long *arem_ptr   = &arem  ; // RE pulse mass
  long *nrem_ptr   = &nrem  ; // RE pulse mass
  long *arew_ptr   = &arew  ; // RE pulse width
  long *nrew_ptr   = &nrew  ; // RE pulse width
  long *arevm_ptr  = &arevm ; // RE pulse mass SD
  long *nrevm_ptr  = &nrevm ; // RE pulse mass SD
  long *arevw_ptr  = &arevw ; // RE pulse width SD
  long *nrevw_ptr  = &nrevw ; // RE pulse width SD
  long *afem_ptr   = &afem  ; // FE pulse mass
  long *nfem_ptr   = &nfem  ; // FE pulse mass
  long *afew_ptr   = &afew  ; // FE pulse width
  long *nfew_ptr   = &nfew  ; // FE pulse width
  long *aetam_ptr  = &aetam ; // pulse mass eta 
  long *aetaw_ptr  = &aetaw ; // pulse width eta
  long *netam_ptr  = &netam ; // pulse mass eta 
  long *netaw_ptr  = &netaw ; // pulse width eta

  // Allocate memory for likelihood 
  likeli = (double *)calloc(1, sizeof(double));

  // Save proposal variances for passing to functions
  sdt   = propsd[6];
  sdrem = propsd[4]; 
  sdrew = propsd[5];
  sdmv  = propsd[2];
  sdwv  = propsd[3];
  sdfem = propsd[7];   // sd of the proposed fixed effect mass distr
  sdfew = propsd[8];   // sd of the proposed fixed effect width distr
  sdetam = propsd[9];  // sd of the proposed eta mass distribution
  sdetaw = propsd[10]; // sd of the proposed eta width distribution


  //----------------------------------------------
  // Create half-life/baseline variance-covariance matrix and decompose
  //----------------------------------------------
  // Allocate matrix memory
  pmd_var = (double **)calloc(2, sizeof(double *));
  for (i=0; i<2; i++) {
    pmd_var[i] = (double *)calloc(2, sizeof(double));
  }

  // Assign proposal values to matrix and calculate covariance
  pmd_var[0][0] = propsd[0];
  pmd_var[1][1] = propsd[1];
  pmd_var[0][1] = pmd_var[1][0] = -0.90 * sqrt(pmd_var[0][0]) * sqrt(pmd_var[1][1]);
  
  // Cholesky decompose the proposal var-covar matrix for b and hl
  pmd_vch = (double **)calloc(2, sizeof(double *));
  for (i = 0; i < 2; i++) {
    pmd_vch[i] = (double *)calloc(2, sizeof(double));
  }

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      pmd_vch[i][j] = pmd_var[i][j];
    }
  }

  if (!cholesky_decomp(pmd_vch, 2)){
    Rprintf("not PSD matrix A\n");
    Rf_error("0");
  }


  //----------------------------------------------------------------
  // Run MCMC
  //----------------------------------------------------------------
  for (i = 0; i < iter; i++) {

    //------------------------------------------------------
    // Run the birth-death algorithm to select pulses 
    //   selection of prior occurs within function based on 
    //   priors->gamma < -0.001
    //------------------------------------------------------
    //if (i < 50) {
    birth_death(list, ts, parms, N, likeli, i, strauss, priors);
    //}

    // Count number of pulses
    // **NOTE**: could remove this traversing of the linked-list by having
    // birth_death() return the number of estimated pulses from this run -- some
    // speed gains likely
    num_node2 = 0;
    new_node = list->succ;
    while (new_node != NULL) {
      new_node = new_node->succ;
      num_node2++;
    }

    //------------------------------------------------------
    // Run gibbs/MH draws
    //------------------------------------------------------

    // 1) Draw the fixed effects   
    //    (Gibbs sampler)
    draw_fixed_effects(list, priors, parms, sdfem, sdfew, afem_ptr, nfem_ptr,
                       afew_ptr, nfew_ptr); 


    // 2) Draw standard deviation of random effects 
    //    (Metropolis Hastings)
    //    Note: log(sd) with uniform prior was suggested by Gelman, 2006
    draw_re_sd(list, priors, parms, sdmv, sdwv, arevm_ptr, nrevm_ptr,
               arevw_ptr, nrevw_ptr); 

    // 3) Draw (kappa from) gamma for the t-distribution var-covar
    //draw_eta(list, parms);
    draw_eta(list, parms, sdetam, sdetaw, aetam_ptr, aetaw_ptr, netam_ptr,
             netaw_ptr);

    // 3) Draw the random effects 
    //    (Metropolis Hastings)
    draw_random_effects(ts, list, parms, N, likeli, sdrem, sdrew, arem_ptr,
                        nrem_ptr, arew_ptr, nrew_ptr); 

    // 4) Draw the pulse locations 
    //    (Metropolis Hastings)
    if (strauss == 1) {
      mh_time_strauss(list, parms, ts, likeli, N, sdt, priors, atime_ptr,
                      ntime_ptr);
    } else {
      mh_time_os(list, parms, ts, likeli, N, sdt, atime_ptr, ntime_ptr,
                 priors->orderstat); 
    }

    // 5) Draw baseline and halflife
    //    (Metropolis-Hastings)
    mh_mu_delta(list, parms, priors, ts, likeli, N, num_node2, pmd_vch,
                adelta_ptr, ndelta_ptr);

    // 6) Draw the model error variance from the inverse Gamma distribution 
    //    (Gibbs) 
    //    Modeling variance NOT precision; precision would be a gamma
    //    distribution via Ken's derivation.  Looked at Week7 of Ed's notes, but
    //    didn't find a clear answer.
    ssq           = error_squared(ts, list, parms, N);
    parms->sigma  = 1 / Rf_rgamma(priors->err_alpha + N / 2, priors->err_beta + 0.5 * ssq);
    parms->lsigma = log(parms->sigma);

    //------------------------------------------------------
    // End MCMC steps
    //------------------------------------------------------


    //------------------------------------------------------
    // Save/Print results
    // TODO: Handle running very few chains or set minimum -- currently if
    // iters < thin, code fails on "attempt to set index 0/0 in SET_VECTOR_ELT"
    //------------------------------------------------------
    // Save to files common/parms
    if ((i >= burnin) && !(i % nthin)) {
      num_node = 0;
      new_node = list->succ;

      // Define output index
      out_i = (i-burnin) / nthin;

      // create matrix for pulse-specific parms in this iteration
      SEXP this_pulse_chain;
      Rf_protect(this_pulse_chain = allocMatrix(REALSXP, num_node2, 8));
      //memset(REAL(this_pulse_chain), 0, num_node*6*sizeof(double)); //not sure how to use this with SEXP matrices
      SEXP dim, dimnames;
      // create 'dim' attribute (dimension of matrix)
      Rf_protect(dim = Rf_allocVector(INTSXP, 2));
      INTEGER(dim)[0] = num_node2; 
      INTEGER(dim)[1] = 8;
      Rf_setAttrib(this_pulse_chain, R_DimSymbol, dim);

      while (new_node != NULL) {

        REAL(this_pulse_chain)[num_node + num_node2*0] = (i + 1);
        REAL(this_pulse_chain)[num_node + num_node2*1] = num_node2;
        REAL(this_pulse_chain)[num_node + num_node2*2] = num_node;
        REAL(this_pulse_chain)[num_node + num_node2*3] = new_node->time;
        REAL(this_pulse_chain)[num_node + num_node2*4] = new_node->theta[0];
        REAL(this_pulse_chain)[num_node + num_node2*5] = new_node->theta[1];
        REAL(this_pulse_chain)[num_node + num_node2*6] = new_node->eta[0];
        REAL(this_pulse_chain)[num_node + num_node2*7] = new_node->eta[1];

        num_node++;
        new_node = new_node->succ;
        //j++;
      }

      // Column names for this_pulse obj in pulse chain --------
      // create names 
      SEXP pulse_chain_names;
      Rf_protect(pulse_chain_names = Rf_allocVector(STRSXP, 8));
      SET_STRING_ELT(pulse_chain_names, 0, Rf_mkChar("iteration"));
      SET_STRING_ELT(pulse_chain_names, 1, Rf_mkChar("total_num_pulses"));
      SET_STRING_ELT(pulse_chain_names, 2, Rf_mkChar("pulse_num"));
      SET_STRING_ELT(pulse_chain_names, 3, Rf_mkChar("location"));
      SET_STRING_ELT(pulse_chain_names, 4, Rf_mkChar("mass"));
      SET_STRING_ELT(pulse_chain_names, 5, Rf_mkChar("width"));
      SET_STRING_ELT(pulse_chain_names, 6, Rf_mkChar("eta_mass"));
      SET_STRING_ELT(pulse_chain_names, 7, Rf_mkChar("eta_width"));
      // assign names to vector
      Rf_protect(dimnames = Rf_allocVector(VECSXP, 2));
      // assign names vector to columns (1) names attribute, leaving rows (0) NULL;
      SET_VECTOR_ELT(dimnames, 1, pulse_chain_names);
      // assign names to list
      Rf_setAttrib(this_pulse_chain, R_DimNamesSymbol, dimnames);

      // Insert 'this_pulse_chain' into list pulse_chain
      SET_VECTOR_ELT(pulse_chains, out_i, this_pulse_chain);

      // Save common parms from iteration (i-burnin)/nthin to SEXP matrix obj -------
      REAL(common)[(out_i) + (out_len)*0] = (i) + 1;
      REAL(common)[(out_i) + (out_len)*1] = num_node2;
      REAL(common)[(out_i) + (out_len)*2] = parms->md[0];
      REAL(common)[(out_i) + (out_len)*3] = parms->theta[0];
      REAL(common)[(out_i) + (out_len)*4] = parms->theta[1];
      REAL(common)[(out_i) + (out_len)*5] = parms->md[1];
      REAL(common)[(out_i) + (out_len)*6] = parms->sigma;
      REAL(common)[(out_i) + (out_len)*7] = parms->re_sd[0];
      REAL(common)[(out_i) + (out_len)*8] = parms->re_sd[1];

      Rf_unprotect(4);

    }


    // Print to STDOUT
    if (verbose & !(i % NNN)) {
      Rprintf("\n\n");
      Rprintf("iter = %d likelihood = %lf\n", i, *likeli);
      Rprintf("mu %.2lf A %.2lf s %.2lf d %.4lf  v %.4le\n", 
              parms->md[0], parms->theta[0], parms->theta[1], parms->md[1],
              parms->sigma);
      Rprintf("pmdvar00 %.2lf pmdvar11 %.2lf pmdvar01 %.2lf \n",
              pmd_var[0][0], pmd_var[1][1], pmd_var[0][1]);
      print_list(list);
      Rprintf("pct rem = %.2lf pct rew = %.2lf pct time = %.2lf\n",
              (double)arem  / (double)nrem, 
              (double)arew  / (double)nrew,
              //(double)afem  / (double)nfem, 
              //(double)afew  / (double)nfew,
              (double)atime / (double)ntime);
      Rprintf("pct md = %.2lf revm = %.2lf revw = %.2lf\n", 
             (double)adelta / (double)ndelta, 
             (double)arevm  / (double)nrevm,
             (double)arevw  / (double)nrevw);
      Rprintf("pct etam = %.2lf etaw = %.2lf\n", 
             (double)aetam / (double)netam, 
             (double)aetaw  / (double)netaw);
    
    }

    //----------------------------------------
    // Check for R user interrupt
    //----------------------------------------
    R_CheckUserInterrupt();


    //------------------------------------------------------
    // Adjust acceptance MH ratios
    //------------------------------------------------------
    if (!(i % 500) && i < 25000 && i > 0) {

      adjust2_acceptance((double) *adelta_ptr / (double) *ndelta_ptr, pmd_var, -0.90);
      adjust_acceptance( (double) *atime_ptr  / (double) *ntime_ptr,  &sdt);
      adjust_acceptance( (double) *arem_ptr   / (double) *nrem_ptr,   &sdrem);
      adjust_acceptance( (double) *arew_ptr   / (double) *nrew_ptr,   &sdrew);
      adjust_acceptance( (double) *arevm_ptr  / (double) *nrevm_ptr,  &sdmv);
      adjust_acceptance( (double) *arevw_ptr  / (double) *nrevw_ptr,  &sdwv);
      adjust_acceptance( (double) *afem_ptr   / (double) *nfem_ptr,   &sdfem);
      adjust_acceptance( (double) *afew_ptr   / (double) *nfew_ptr,   &sdfew);
			adjust_acceptance( (double) *aetam_ptr  / (double) *netam_ptr,  &sdetam);
			adjust_acceptance( (double) *aetaw_ptr  / (double) *netaw_ptr,  &sdetaw);

      adelta = ndelta = 0;
      atime  = ntime  = 0;
      arem   = nrem   = 0;
      arew   = nrew   = 0;
      arevm  = nrevm  = 0;
      arevw  = nrevw  = 0;
      afem   = nfem   = 0;
      afew   = nfew   = 0;
      aetam  = netam  = 0;
      aetaw  = netaw  = 0;

      // Cholesky decompose the three proposal var-covar matrices
      for (k = 0; k < 2; k++) {
        for (l = 0; l < 2; l++) {
          pmd_vch[k][l] = pmd_var[k][l];
        }
      }

      if (!cholesky_decomp(pmd_vch, 2)){
        Rprintf("pmd not PSD matrix\n");
        Rf_error("0");
      }

    } 

  } // End of MCMC loop 

  // Close files and free allocated memory
  //fclose(common);
  //fclose(parm);
  free(likeli);
  PutRNGstate();

} 




//-----------------------------------------------------------------------------
//
// mh_time_strauss: 
//   this runs the M-H draws for each individual pulse location using a
//   Strauss process prior - can be set to a hard core by setting gamma = 0
//
//   NOTE: Estimated on natural scale of location/time
//
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that
//                           exist
//     Common_parms *parms - the current values of the common
//                           parameters
//     double **ts         - this is the matrix of observed data (a
//                           column of times and a column of log(concentration)
//     double *like        - the current value of the likelihood
//     int N               - the number of observations in **ts
//     double v            - the proposal variance for pulse location
//     Priors *priors      - For priors->gamma, the repulsion parameter for
//                           Strauss process (hard core: gamma=0); and 
//                           priors->range, the range of repulsion (R) for
//                           Strauss process 
//
//   RETURNS: 
//     None                - all updates are made internally
// 
//-----------------------------------------------------------------------------
void mh_time_strauss(Node_type *list, 
                     Common_parms *parms, 
                     double **ts, 
                     double *like, 
                     int N, 
                     double sd, 
                     Priors *priors,
                     long *atime,
                     long *ntime) {

  // Declarations ----------------------
  int i;                    // Generic counter
  int sum_s_r_proposal;     // sum[S(R)] using proposed location
  int sum_s_r_current;      // sum[S(R)] using current location
  double alpha;             // Acceptance ratio
  double plikelihood;       // Likelihood using proposed value
  double ptime;             // Proposal location
  double like_ratio;        // Likelihood portion of acceptance ratio
  double prior_ratio;       // Prior portion of acceptance ratio
  double l_rho;             // Log of rho (acceptance ratio prior to min() portion)
  double current_time;      // current location
  double *tmp_mean_contrib; // Temp array of mean contribution for each pulse
  Node_type *node;          // Pointer to node in linked list


  // Allocate Memory 
  tmp_mean_contrib = (double *)calloc(N, sizeof(double));

  // Set pointer to first node in linked list 
  node = list->succ;

  while (node != NULL) {

    // Increase denominator of acceptance rate for time 
    (*ntime)++;

    // Compute proposal time 
    ptime = Rf_rnorm(node->time, sd); 

    // Calculate sum_s_r for proposal value and current value
    sum_s_r_proposal = calc_sr_strauss(ptime, list, node, priors);
    sum_s_r_current  = calc_sr_strauss(node->time, list, node, priors);

    // if valid proposal time, run alpha and accept/reject code
    if (ptime <= fitend && ptime > fitstart) { 

      // Calculate prior ratio - 
      //      gamma^(sum(s_r*)_proposal - sum(s_r)_current)
      prior_ratio = pow(priors->gamma, sum_s_r_proposal - sum_s_r_current);

      // if prior ratio is 0 (~EPS), set l_rho to 0, 
      // else calculate it 
      // (note: necessary because log(0) = NaN, like_ratio on log scale)
      if (fabs(prior_ratio) < EPS) {

        l_rho = -1e300;

      } else {

        // Save current time and set its time to proposed value
        current_time = node->time;
        node->time   = ptime;

        // Save mean_contrib of this pulse 
        for (i = 0; i < N; i++) {
          tmp_mean_contrib[i] = node->mean_contrib[i];
        }

        // Recalculate mean_contrib of this pulse assuming proposed value 
        mean_contribution(node, ts, parms, N);

        // Calculate the likelihood under the proposed value 
        // returns log-likelihood 
        plikelihood = likelihood(list, ts, parms, N, list);

        // Calculate the likelihood ratio 
        like_ratio = plikelihood - *like;
        l_rho      = (log(prior_ratio) + like_ratio);

      }

      // Calculate log rho; set alpha equal to min(0, log rho) 
      // x ? y:z is equivalent to R's ifelse, y = true, z = false 
      alpha = (0 < l_rho) ? 0 : l_rho;

      // Accept/Reject
      if (log(Rf_runif(0, 1)) < alpha) {

        // If log U < log rho, we accept proposed value. Increase
        // acceptance count by one and set likelihood equal to
        // likelihood under proposal 
        (*atime)++;
        *like = plikelihood;

      } else if (fabs(prior_ratio) > EPS) {

        // Otherwise, we reject, and if prior_ratio was != 0, we have to
        // revert back to current values. if prior_ratio ~= 0, we didn't
        // change the values -- just set l_rho

        // Set pulse's time back to current time 
        node->time = current_time;

        // Set mean_contrib of this pulse back to current value 
        for (i = 0; i < N; i++) {
          node->mean_contrib[i] = tmp_mean_contrib[i];
        }

      }

    } // End of if time is feasible statement

    // Advance one pulse
    node = node->succ;

  } // End of loop through pulses

  // Free Memory
  free(tmp_mean_contrib);

}




//-----------------------------------------------------------------------------
//
// mh_time_os: 
//   this runs the M-H draws for each individual pulse location
//
//   NOTE: Estimated on natural scale of location/time
// 
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that
//                           exist
//     Common_parms *parms - the current values of the common
//                           parameters
//     double **ts         - this is the matrix of observed data (a
//                           column of times and a column of
//                           log(concentration)
//     double *like        - the current value of the likelihood
//     int N               - the number of observations in **ts
//     double v            - the proposal variance for pulse location
//     int mmm             - Order stat to use for prior on location
//                           (priors->orderstat)
//
//   RETURNS: 
//       None                - all updates are made internally
// 
//-----------------------------------------------------------------------------
void mh_time_os(Node_type *list, 
                Common_parms *parms, 
                double **ts, 
                double *like,
                int N, 
                double sd,
                long *atime,
                long *ntime,
                int mmm) {

  int i;                     // Generic counter
  Node_type *node;           // Pointer for current list of pulses
  double alpha;              // min(0, log(rho)) - acceptance ratio
  double plikelihood;        // Likelihood using proposal location
  double ptime;              // Proposed location
  double like_ratio;         // Likelihood portion of acceptance ratio
  double prior_ratio;        // Prior portion of acceptance ratio
  double rho;                // acceptance ratio prior to calc'ing alpha
  double current_time;       // Current location
  double *curr_mean_contrib; // Mean contrib w/ current location
  double time_diff1;         // For curr location and prev pulse
  double time_diff1_new;     // For proposed location and prev pulse
  double time_diff2;         // For curr location and succ pulse
  double time_diff2_new;     // For proposed location and succ pulse

  // Allocate Memory 
  curr_mean_contrib = (double *)calloc(N, sizeof(double));

  node = list->succ;

  while (node != NULL) {

    // Increase denominator of acceptance rate for time 
    (*ntime)++;

    // Compute proposal time 
    ptime = Rf_rnorm(node->time, sd);

    // Only proceed if our proposed time is reasonable 
    if (ptime <= fitend && ptime > fitstart) {

      if (node->pred != NULL) {
        // If we're at the first pulse, use these as the first
        // numerator/denominator of the prior ratio 
        time_diff1_new = ptime - node->pred->time;
        time_diff1 = node->time - node->pred->time;
      } else {
        // Otherwise, use these as the first numerator/denominator of
        // the prior ratio 
        time_diff1_new = ptime - fitstart;
        time_diff1 = node->time - fitstart;
      }

      if (node->succ != NULL) {
        // If we're at the last pulse, use these as the second
        // numerator/ denominator of the prior ratio 
        time_diff2_new = node->succ->time - ptime;
        time_diff2 = node->succ->time - node->time;
      } else {
        // Otherwise, use these as the second numerator/denominator 
        time_diff2_new = fitend - ptime;
        time_diff2_new = fitend - node->time;
      }

      // Combine it all for the prior ratio 
      prior_ratio = (mmm-1) * 
                    log((time_diff1_new * time_diff2_new) / 
                        (time_diff1 * time_diff2));

      // Save current time 
      current_time = node->time;

      // Set pulse's time to the proposal value 
      node->time = ptime;

      // Save mean_contrib of this pulse 
      for (i = 0; i < N; i++) {
        curr_mean_contrib[i] = node->mean_contrib[i];
      }

      // Recalculate mean_contrib of this pulse assuming proposed value 
      mean_contribution(node, ts, parms, N);

      // Calculate the likelihood under the proposed value returns
      // log-likelihood 
      plikelihood = likelihood(list, ts, parms, N, list);

      // Calculate the likelihood ratio 
      like_ratio = plikelihood - *like;

      // Calculate log rho; set alpha equal to min(0,log rho) 
      // x ? y:z is equivalent to R's ifelse, y = true, z = false 
      alpha = (0 < (rho = (prior_ratio + like_ratio))) ? 0 : rho;

      if (log(Rf_runif(0, 1)) < alpha) {
        // If log U < log rho, we accept proposed value. Increase acceptance
        // count by one and set likelihood equal to likelihood under proposal 
        (*atime)++;
        *like = plikelihood;
      } else {
        // Otherwise, we reject, and we have to revert back to current
        // values 

        // Set pulse's time back to current time 
        node->time = current_time;

        // Set mean_contrib of this pulse back to current value 
        for (i = 0; i < N; i++) {
          node->mean_contrib[i] = curr_mean_contrib[i];
        }
      }

    } 
    // End of if time is feasible statement

    // Advance one pulse
    node = node->succ;

  }
  // End of loop through pulses

  // Free Memory
  free(curr_mean_contrib);

}




//-----------------------------------------------------------------------------
//
// mh_mu_delta: 
//   this runs the M-H draw for baseline and halflife (drawn together)
//
//   NOTE: Estimated on natural scale of baseline and half-life
// 
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that
//                           exist
//     Common_parms *parms - the current values of the common
//                           parameters
//     Priors *priors      - the parameters of the prior distributions
//     double **ts         - this is the matrix of observed data (a
//                           column of times and a column of
//                           log(concentration)
//     double *like        - the current value of the likelihood
//     int N               - the number of observations in **ts
//     int num_node        - the current number of pulses
//     double **var        - the proposal variance-covariance matrix
//                           for baseline and halflife
//
//   RETURNS: 
//     None                - all updates are made internally
//
//-----------------------------------------------------------------------------
void mh_mu_delta(Node_type *list, 
                 Common_parms *parms, 
                 Priors *priors, 
                 double **ts, 
                 double *like, 
                 int N, 
                 int num_node, 
                 double **var,
                 long *adelta,
                 long *ndelta) {

  int i;                // Generic counter
  int j;                // Generic counter
  int k;                // Generic counter
  Node_type *node;      // Point to pulses linklist
  double alpha;         // min(0,log(rho)), acceptance ratio
  double plikelihood;   // likelihood with proposed value
  double *pmd;          // array of proposed values (bl/hl)
  double like_ratio;    // Diff btwn curr and proposed likelihoods
  double priorb_old;    // Portion of accept ratio due to curr value/prior (bl)
  double priorb_new;    // Portion of accept ratio due to prev value/prior (bl)
  double priorh_old;    // Portion of accept ratio due to curr value/prior (hl)
  double priorh_new;    // Portion of accept ratio due to prev value/prior (hl)
  double prior_ratio;   // Diff btwn curr and proposed priors
  double currentmd[2];  // Array of current values of bl/hl
  double **currentmc;   // Current mean contrib for each pulse
  double logrho;        // Acceptance ratio prior to calc'ing alpha
  double current_decay; // Current decay rate

  // Allocate memory 
  currentmc = (double **)calloc(num_node, sizeof(double *));
  for (i = 0; i < num_node; i++) {
    currentmc[i] = (double *)calloc(N, sizeof(double));
  }
  pmd = (double *)calloc(2, sizeof(double));

  // Increase denominator of acceptance rate for b and hl 
  (*ndelta)++;

  // Draw proposal values for b and hl 
  rmvnorm(pmd, var, 2, parms->md, 1); 

  // Only proceed if we draw reasonable values 
  if (pmd[0] > 0 && pmd[1] > 3) {

    // Compute ratio of prior densities 
    priorb_old  = parms->md[0] - priors->meanbh[0];
    priorb_old *= 0.5 * priorb_old / priors->varbh[0];
    priorb_new  = pmd[0] - priors->meanbh[0];
    priorb_new *= 0.5 * priorb_new / priors->varbh[0];
    priorh_old  = parms->md[1] - priors->meanbh[1];
    priorh_old *= 0.5 * priorh_old / priors->varbh[1];
    priorh_new  = pmd[1] - priors->meanbh[1];
    priorh_new *= 0.5 * priorh_new / priors->varbh[1];

    prior_ratio = priorb_old + priorh_old - priorb_new - priorh_new;

    // Save current values of b and hl and set new current values equal to
    // proposal values 
    for (k = 0; k < 2; k++) {
      currentmd[k] = parms->md[k];
      parms->md[k] = pmd[k];
    }

    // Save current decay rate; calculate new current decay rate based on
    // proposal value of hl 
    current_decay = parms->decay;
    parms->decay  = log(2) / parms->md[1];

    // Save current mean_contrib for each pulse; calculate new current
    // mean_contrib for each pulse based on proposed values of b and hl */
    i = 0;
    node = list->succ;
    while (node != NULL) {
      for (j = 0; j < N; j++) {
        currentmc[i][j] = node->mean_contrib[j];
      }
      mean_contribution(node, ts, parms, N);
      i++;
      node = node->succ;
    }

    // Calculate proposed likelihood and then calculate likelihood ratio 
    plikelihood = likelihood(list, ts, parms, N, list);
    like_ratio = plikelihood - *like;

    // Calculate log rho; set alpha equal to min(0, log rho) 
    alpha = (0 < (logrho = (prior_ratio + like_ratio))) ? 0:logrho;

    if (log(Rf_runif(0, 1)) < alpha) {
      // If log U < log rho, increase acceptance rate by 1 and set the
      // likelihood equal to the likelihood under proposed values 
      (*adelta)++;
      *like = plikelihood;
    } else {
      // Otherwise, we need to revert values back to their current state 

      // Set b and hl back equal to current values
      parms->md[0] = currentmd[0];
      parms->md[1] = currentmd[1];

      // Set mean_contrib back equal to current values for each pulse
      i = 0;
      node = list->succ;
      while (node != NULL) {
        for (j = 0; j < N; j++) {
          node->mean_contrib[j] = currentmc[i][j];
        }
        i++;
        node = node->succ;
      }

      // Set decay rate back equal to current value
      parms->decay = current_decay;
    } 
    // End of if else statement

  }

  // Free memory 
  for (i = 0; i < num_node; i++) {
    free(currentmc[i]);
  }
  free(currentmc);
  free(pmd);

}




//-----------------------------------------------------------------------------
// 
// draw_fixed_effects (OLD - lognormal): 
//   this runs the M-H draw for overall mean mass and width
// 
//   NOTE: Estimated on log scale, relative to natural scales of mass and
//     width.
//
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that
//                           exist      
//     Priors *priors      - the parameters of the prior distributions          
//     Common_parms *parms - the current values of the common
//                           parameters        
//     double v1           - the proposal variance for overall mean
//                           mass        
//     double v2           - the proposal variance for overall mean
//                           width       
// 
//   RETURNS: 
//     None                - all updates are made internally
// 
//-----------------------------------------------------------------------------
//void draw_fixed_effects(Node_type *list, 
//                        Priors *priors, 
//                        Common_parms *parms
//                        ) {
//
//  // Declarations
//  int j;           // Generic counter
//  int numnode;     // Count of pulses
//  double gmean;    // Mean for Gibbs sampler 
//  double gvar;     // Var for Gibbs sampler
//  double resum;    // Sum of pulse individual masses/widths
//  Node_type *node; // Pointer to linklist of pulses
//
//
//  // Draw new values of mu_a and mu_w via Gibbs Sampler 
//  for (j = 0; j < 2; j++) {
//
//    resum   = 0.0;
//    numnode = 0;
//    node    = list->succ;
//
//    // This while loop counts the pulses and gets the sum of r.e. 
//    while (node != NULL) {
//      numnode++;
//      resum += log(node->theta[j]);
//      node   = node->succ;
//    } 
//
//    // I added this on 11 Sep 2012; numnode is 0 to start with and this
//    // causes gvar to be very large, sometimes leading to values of mua and
//    // muw that are too large or too small; this in turn leads to later
//    // extreme values of A-ik and sigma-p-ik which causes precision issues
//    // in mean-contrib 
//    if (numnode == 0) {
//      numnode = 1;
//    }
//
//    gmean  = (priors->fe_variance[j] * resum);
//    gmean += (parms->re_sd[j] * parms->re_sd[j] * priors->fe_mean[j]);
//    gmean /= (priors->fe_variance[j] * numnode) + 
//             (parms->re_sd[j] * parms->re_sd[j]);
//    gvar   = priors->fe_variance[j] * parms->re_sd[j] * parms->re_sd[j];
//    gvar  /= (priors->fe_variance[j] * numnode) + 
//             (parms->re_sd[j] * parms->re_sd[j]);
//
//    parms->theta[j] = Rf_rnorm(gmean, sqrt(gvar));
//  }
//
//}



//-----------------------------------------------------------------------------
// START OF draw_fixed_effects SUBROUTINE*/
//
//
// draw_fixed_effects: this runs the Gibbs sampler draw for subject-specific
// mean masses and widths;
// ARGUMENTS: Subject_type *sublist; this is the current list of subjects;
// Priors *priors; the current values of the prior parameters;
// Common_parms *parms; the current values of the common parameters;
// unsigned long *seed; seed values needed for the randon number generator;
// RETURNS: None; all updates are made internally
//
//********************************************************************/
//
// VARIABLE DEFINITIONS
// j: generic counter
// numnode: counter of pulses
//
// gvar: variance of distribution used in Gibbs sampler
// *subject: current list of subjects and their characteristics
// *subnode: for a given subject, list of pulses and their characteristics
//
// SUBROUTINES USED
// kiss: found in randgen.c; draws from U(0,1) distribution
// rnorm: found in randgen.c; draws from the normal distribution
//------------------------------------------------------------------------------
void draw_fixed_effects(Node_type *list, 
                        Priors *priors, 
                        Common_parms *parms, 
                        double sdfem, 
                        double sdfew,
                        long *afem, long *nfem, long *afew, long *nfew ) {  
  int j, newint, oldint, numnode;
  long *accept_counter;
  double normalizing_ratio, acceptance_ratio, theta[2], old_prior, new_prior,
         prior_ratio, alpha, psum_old,
         psum_new, prop_ratio, stdxold, stdxnew;
  Node_type *node;

  accept_counter = (long *)calloc(2, sizeof(long));

  // Acceptance counters
  accept_counter[0]  = *afem;
  accept_counter[1]  = *afew;
  // Iteration counters
  (*nfem)++;
  (*nfew)++;

  // Proposed values
  theta[0] = Rf_rnorm(parms->theta[0], sdfem);
  theta[1] = Rf_rnorm(parms->theta[1], sdfew);

  // Draw a new pair of pulse mass and width on the subject level
  for (j = 0; j < 2; j++) {

    if (theta[j]>0) {

      // Prior Ratio  -- DONE
      old_prior    = (parms->theta[j] - priors->fe_mean[j]) * 
                     (parms->theta[j] - priors->fe_mean[j]); 
      new_prior    = (theta[j] - priors->fe_mean[j]) * 
                     (theta[j] - priors->fe_mean[j]); 
      prior_ratio  = (old_prior - new_prior) / (2 * priors->fe_variance[j]);

      // likelihood ratio (Proposal Ratio)
      psum_old = 0;
      psum_new = 0;
      newint   = 0;
      oldint   = 0;

      numnode = 0;
      node = list->succ;

      // 'likelihood' ratio -- Ratio of p(alpha|mu, nu, kappa)
      while (node != NULL) {
        psum_old += (node->theta[j] - parms->theta[j]) * 
                    (node->theta[j] - parms->theta[j]) * node->eta[j];
        psum_new += (node->theta[j] - theta[j]) * 
                    (node->theta[j] - theta[j]) * node->eta[j];

        // Normalizing constants
        stdxnew   = theta[j]        * sqrt(node->eta[j]) / parms->re_sd[j];
        stdxold   = parms->theta[j] * sqrt(node->eta[j]) / parms->re_sd[j];
        oldint += Rf_pnorm5(stdxold, 0, 1, 1.0, 1.0); // second 1.0 does the log xform for us 
        newint += Rf_pnorm5(stdxnew, 0, 1, 1.0, 1.0); // first 1.0 says to use lower tail

        node = node->succ;
        numnode++;
      }

      prop_ratio = 0.5 / (parms->re_sd[j] * parms->re_sd[j]) * 
                   (psum_old - psum_new);
      normalizing_ratio = oldint - newint; 

      acceptance_ratio = prior_ratio + prop_ratio + normalizing_ratio;
      alpha = (0 < acceptance_ratio) ? 0 : acceptance_ratio;

      // If log(U) < log rho, accept the proposed value
      // Increase acceptance count by 1
      if (log(Rf_runif(0, 1)) < alpha) {
        accept_counter[j]++;
        parms->theta[j] = theta[j];
      } 

    } 

  } // end of loop through mass & width 

  *afem = accept_counter[0];
  *afew = accept_counter[1];
  free(accept_counter);

}





//-----------------------------------------------------------------------------
//
// draw_re_sd: 
//   this runs the M-H draw for overall standard deviations of mass and
//   width
// 
//   NOTE: This function was updated due to a change in priors on alpha_i and
//         omega_i to a truncated t.  This re-introduced an occasional
//         sticky-point issue where the re_sd (esp. for width) goes to 0 and
//         stays there for 1000 to 100000 iterations.  Changing Unif(0,max) to
//         Unif(0.1,max) seems to alleviate the issue, but doesn't fix the core
//         (unidentified) problem. Is there a numerical issue?  Perhaps in the
//         normalizing constants? 
//
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that exist      
//     Priors *priors      - the parameters of the prior distributions          
//     Common_parms *parms - the current values of the common parameters        
//     double v1           - the proposal variance for overall st-dev of mass   
//     double v2           - the proposal variance for overall st-dev of width  
//
//   RETURNS: 
//     None                - all updates are made internally
//
//-----------------------------------------------------------------------------
void draw_re_sd(Node_type *list,
                Priors *priors, 
                Common_parms *parms, 
                double v1, 
                double v2,
                long *arevm, long *nrevm, long *arevw, long *nrevw) {
  int j, 
      num_pulses;
  long *accept_counter;   
  double *new_sd, 
         stdx_old, 
         stdx_new, 
         new_int, 
         old_int, 
         first_part,
         second_part, 
         third_part, 
         log_rho;
  Node_type *node;

  // Allocate Memory
  new_sd         = (double *)calloc(2, sizeof(double));
  accept_counter = (long *)calloc(2, sizeof(long));

  // Add 1 to the counters for acceptance rates of sigma_a and sigma_w
  (*nrevm)++;
  (*nrevw)++;

  // Assign current acceptance counts to temporary vector
  accept_counter[0] = *arevm;
  accept_counter[1] = *arevw;

  // Draw proposed values for sigma_a and sigma_w
  new_sd[0] = Rf_rnorm(parms->re_sd[0], v1);
  new_sd[1] = Rf_rnorm(parms->re_sd[1], v2);

  // Accept or Reject sigma_a, then accept or reject for sigma_w
  for (j = 0; j < 2; j++) {

    // We only can accept the proposed value if it is positive
    if (new_sd[j] > 0 && new_sd[j] < priors->re_sdmax[j]) {

      node = list->succ;
      num_pulses = 0;
      third_part = 0;
      old_int = new_int = 0;
      while (node != NULL) {

        // Normalizing constants
        stdx_old   = parms->theta[j] / ( parms->re_sd[j] / sqrt(node->eta[j]) );
        stdx_new   = parms->theta[j] / ( new_sd[j]       / sqrt(node->eta[j]) );
        new_int   += Rf_pnorm5(stdx_new, 0, 1, 1.0, 1.0);
        old_int   += Rf_pnorm5(stdx_old, 0, 1, 1.0, 1.0);

        // Count pulses
        num_pulses++;
        // 3rd 'part' of acceptance ratio
        third_part += node->eta[j] * 
                      (node->theta[j] - parms->theta[j]) *
                      (node->theta[j] - parms->theta[j]);

        // Next pulse
        node = node->succ;
      }

      // 1st and 2nd 'parts' of acceptance ratio
      first_part  = (num_pulses) * (log(parms->re_sd[j]) - log(new_sd[j]));
      second_part = 0.5 * 
        ((1 / (parms->re_sd[j] * parms->re_sd[j])) - (1 / (new_sd[j] * new_sd[j])));

      // Compute log rho, and set alpha equal to min(log rho,0)
      log_rho = old_int - new_int + first_part + second_part * third_part;
      log_rho = fmin(0, log_rho);

      // If log(U) < log rho, accept the proposed value
      if (log(Rf_runif(0, 1)) < log_rho) {
        accept_counter[j]++;
        parms->re_sd[j] = new_sd[j];
      }

    }

  } 

  // Set acceptance count equal to temp vector components
  *arevm = accept_counter[0];
  *arevw = accept_counter[1];

  free(new_sd);
  free(accept_counter);

}



//-----------------------------------------------------------------------------
// 
// draw_random_effects: 
//   this runs the M-H draw for individual pulse masses and widths
//  
//   NOTE: now uses a truncated t prior on alpha_i and omega_i
// 
//   ARGUMENTS: 
//     double **ts         - this is the matrix of observed data (a
//                           column of times and a column of log(concentration) 
//     Node_type *list     - this is the current list of pulses that exist                                             
//     Common_parms *parms - the current values of the common parameters
//     int N               - the number of observations in **ts
//     double *like        - the current value of the likelihood
//     double v1           - the proposal variance for individual pulse masses                                         
//     double v2           - the proposal variance for individual pulse widths                                         
// 
//   RETURNS: 
//     None                - all updates are made internally
// 
// 
//-----------------------------------------------------------------------------
void draw_random_effects(double **ts, 
                         Node_type *list, 
                         Common_parms *parms, 
                         int N, 
                         double *like,
                         double v1, 
                         double v2, 
                         long *arem, long *nrem, 
                         long *arew, long *nrew) {
  //double current_like;
  int i;                  // Generic counters
  int j;                  // Generic counters
  double logrho;          // Acceptance ratio prior to min(0,logrho)
  double prior_old;       // Prior portion of accept. ratio due to current value
  double prior_new;       // Prior portion of accept. ratio due to proposed value
  double old_val;         // For holding old value of mass/width
  long *accept_counter; // Internal array for counting acceptances
  double *pRE;            // Array of proposed RE mass/width value
  double prior_ratio;     // Portion of accept. ratio due to priors
  double like_ratio;      // Portion of accept. ratio due to likelihoods
  double alpha;           // Acceptance ratio after min(0,logrho)
  double plikelihood;     // Likelihood with proposed value
  double *old_contrib;    // For holding mean contrib using current value
  Node_type *node;        // Pointer to linklist of pulses


  // Go to start of pulses 
  node = list->succ;

	// Set acceptance counts equal to temporary vector
	accept_counter    = (long *)calloc(2, sizeof(long));
	accept_counter[0] = *arem;
	accept_counter[1] = *arew;

  // Go through each existing pulse 
  while (node != NULL) {

    // Allocate Memory
    old_contrib = (double *)calloc(N, sizeof(double));
    pRE         = (double *)calloc(2, sizeof(double));

    // Increase the denominators of the acceptance rates
    (*nrem)++;
    (*nrew)++;

    // Draw proposed values of current pulse's mass and width
    pRE[0] = Rf_rnorm(node->theta[0], v1);
    pRE[1] = Rf_rnorm(node->theta[1], v2);

    // Determine if we accept or reject proposed pulse mass then determine
    // if we accept or reject proposed pulse width
    for (j = 0; j < 2; j++) {
      //if ((j == 0 && pRE[0] > 0) | (j == 1 && pRE[1] > 5)) {
      
      if (pRE[j] > 0) {

        // Compute the log of the ratio of the priors
        prior_old = node->theta[j] - parms->theta[j];
        prior_old *= 0.5*prior_old;
        prior_new = pRE[j] - parms->theta[j];
        prior_new *= 0.5*prior_new;
        prior_ratio = prior_old - prior_new;
        prior_ratio /= parms->re_sd[j];
        prior_ratio /= parms->re_sd[j];

        // Save the current value of mass/width
        old_val = node->theta[j];

        // Set the pulse's mass/width equal to the proposed value
        node->theta[j] = pRE[j];

        // Save the mean_contrib for that pulse
        for (i = 0; i < N; i++) {
          old_contrib[i] = node->mean_contrib[i];
        }

        // Recalculate that pulse's mean_contrib assuming proposed mass/width 
        mean_contribution(node, ts, parms, N);

        // Calculate likelihood assuming proposed mass/width 
        plikelihood = likelihood(list, ts, parms, N, list);

        like_ratio = plikelihood - *like;

        // Compute the log of the ratio between the two likelihoods

        // Calculate log rho; set alpha equal to min(0,log rho) 
        alpha = (0 < (logrho = (prior_ratio + like_ratio))) ? 0:logrho;

        // If log U < log rho, accept the proposed value, increase acceptance
        // counter 
        if (log(Rf_runif(0, 1)) < alpha) {
          accept_counter[j]++;
        } else {
        // Otherwise, reject the proposed value, set pulse's mass/width back to 
        // saved current value, and set pulse's mean_contrib back to saved value

          node->theta[j] = old_val;
          for (i = 0; i < N; i++) {
            node->mean_contrib[i] = old_contrib[i];
          }

        }
      } // end of loop through mass & width 
    }

    // Advance to next pulse
    node = node->succ;
    free(pRE);
    free(old_contrib);

  } // end of loop through pulses

  // Set counter equal to temporary vector values
  *arem = accept_counter[0];
  *arew = accept_counter[1];
  // free memory
  free(accept_counter);

}



//-----------------------------------------------------------------------------
// 
// draw_eta: 
//   This subroutine draws a kappa for the scaling of the SD in the truncated-t
//   prior on alpha_i/omega_i (individual mass/width). 
// 
//-----------------------------------------------------------------------------
void draw_eta(Node_type *list, 
              Common_parms *parms, 
              double sdm, double sdw,
              long *aetam, long *aetaw, long *netam, long *netaw) {
  int j;
  long *accept_counter, num_node; 
  double prior_ratio, re_ratio, old_gamma, new_gamma, stdold, stdnew,
         re_old, re_new, alpha, temp;
  double *proposed_eta;
  Node_type *node;

  proposed_eta      = (double *)calloc(2, sizeof(double));
  accept_counter    = (long *)calloc(2, sizeof(long));
  accept_counter[0] = *aetam;
  accept_counter[1] = *aetaw;

  // first pulse mass and then pulse width
  num_node = 0;
  node = list->succ;
  while (node != NULL) {
    (*netam)++;
    (*netaw)++;
    num_node++;

    // draw the new eta
    proposed_eta[0] = Rf_rnorm(node->eta[0], sdm);
    proposed_eta[1] = Rf_rnorm(node->eta[1], sdw);

    for (j = 0; j < 2; j++) {

      if (proposed_eta[j] > 0) {

        // Shape, scale parameterized: 
        //    https://github.com/mmulvahill/r-source/blob/trunk/src/nmath/dgamma.c
        //    https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Distribution-functions
        old_gamma = Rf_dgamma(node->eta[j], 2, 0.5, 0); 
        new_gamma = Rf_dgamma(proposed_eta[j], 2, 0.5, 0);

        prior_ratio  = log(new_gamma) - log(old_gamma);
        stdold       = (node->theta[j]) / (parms->re_sd[j] / sqrt(node->eta[j]));
        stdnew       = (node->theta[j]) / (parms->re_sd[j] / sqrt(proposed_eta[j]));
        re_old       = node->theta[j] - parms->theta[j];
        re_old      *= 0.5 * re_old * node->eta[j];
        re_new       = node->theta[j] - parms->theta[j];
        re_new      *= 0.5 * re_new * proposed_eta[j];
        re_ratio     = re_old - re_new;
        re_ratio    /= parms->re_sd[j];
        re_ratio    /= parms->re_sd[j];
        re_ratio    += Rf_pnorm5(stdold, 0, 1, 1.0, 1.0) -  // second 1.0 does the log xform for us 
                       Rf_pnorm5(stdnew, 0, 1, 1.0, 1.0) -  // first 1.0 says to use lower tail      
                       0.5 * log(node->eta[j]) + 0.5 * log(proposed_eta[j]); // the 1/2pi term in
                                                                             // normal distirbution
        alpha = (0 < (temp = (prior_ratio + re_ratio))) ? 0 : temp;

        // If log U < log rho, accept the proposed value, increase acceptance
        // counter 
        if (log(Rf_runif(0, 1)) < alpha) {
          accept_counter[j]++;
          node->eta[j] = proposed_eta[j];
        }

      }

      // NOTE: eta set to 1 for debugging
      //node->eta[j] = 1;

    }

    node = node->succ;

  }
  *aetam = accept_counter[0];
  *aetaw = accept_counter[1];
  free(accept_counter);
  free(proposed_eta);
  
}



//-----------------------------------------------------------------------------
// 
// error_squared: 
//   This subroutine calculates the sum of squared error. It uses the list of
//   pulses and the common parameters to evaluate the likelihood and
//   calculate the sum of the squared error between observed concentration
//   and expected under the current parameters
// 
//   ARGUMENTS: 
//     double **ts         - this is the matrix of observed data (a column
//                           of times and a column of log(concentration)
//     Node_type *list     - this is the current list of pulses that exist
//     Common_parms *parms - the current values of the common parameters
//     int N               - the number of observations in **ts
// 
//   RETURNS: 
//     double ssq          - this value is needed for one of the parameters
//                             of the posterior distribution of model error
//                             variance
// 
//-----------------------------------------------------------------------------
double error_squared(double **ts, 
                     Node_type *list, 
                     Common_parms *parms, 
                     int N) {

  // Declarations
  int i;        // Generic counter
  double *mean; // mean concentration w/ current parms and pulses
  double ssq;   // sum of squared differences btwn log(conc) and expected values

  mean = mean_concentration(list, parms, N, list, ts);
  ssq  = 0;
  for (i = 0; i < N; i++) {
    ssq += (ts[i][1] - mean[i]) * (ts[i][1] - mean[i]);
  }

  free(mean);

  return ssq;

}




//-----------------------------------------------------------------------------
//
// adjust_acceptance: 
//   this adjusts the proposal variances based on the inputted acceptance
//   rate and proposal variance. If the acceptance rate is too high or too
//   low, proposal variance will be adjusted.
//
//   ARGUMENTS: 
//     double x  - the inputted acceptance rate; usually inputted as the
//                 acceptance counter divided by the attempt counter
//     double *X - the current proposal variance
//
//   RETURNS: 
//     None      - update to the proposal variance is made internally
//
//-----------------------------------------------------------------------------
void adjust_acceptance(double x, double *X) {

  double y; // new proposal variance based on inputs

  y = 1. + 1000. * (x-.35) * (x-.35) * (x-.35);

  if (y < .9) {
    y = .9;
  }

  if (y > 1.1) {
    y = 1.1;
  }

  *X *= y;

}




//-----------------------------------------------------------------------------
// 
// adjust2_acceptance: 
//   this adjusts the proposal variance-covriance matrix for baseline and
//   halflife based on the inputted acceptance rate and proposal matrix. If
//   the acceptance rate is too high or too low, proposal variance will be
//   adjusted.
// 
//   ARGUMENTS: 
//     double x    - the inputted acceptance rate; usually inputted as the
//                   acceptance counter divided by the attempt counter
//     double **X  - the current proposal variance-covariance matrix
//     double corr - the correlation between the two proposal variances
// 
//   RETURNS: 
//     None        - update to the proposal variance is made internally
// 
// 
//-----------------------------------------------------------------------------
void adjust2_acceptance(double x, double **X, double corr) {

  double y; // new diagonal elements of proposal var-covar matrix

  y = 1. + 1000. * (x - .25) * (x - .25) * (x - .25);

  if (y < .90) {
    y = .90;
    X[0][0] *= y;
    X[1][1] *= y;
    X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
  }

  if (y > 1.1) {
    y        = 1.1;
    X[0][0] *= y;
    X[1][1] *= y;
    X[0][1]  = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
  }

}




//------------------------------------------------------------------------------
// END OF FILE
//------------------------------------------------------------------------------

