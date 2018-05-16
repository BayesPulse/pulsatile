//-----------------------------------------------------------------------------
//
// FILE: birth_death.c 
// AUTHOR: Matt Mulvahill 
//         (and many before me: T Johnson, N Carlson, K Horton, Karen )
//
// DESCRIPTION: 
//   Contains birth-death algorithm, expanded to spatial birth-death with the
//   option to use the older order-statistic prior approach.  To use the order
//   statistic, priors->gamma should be less than some very small negative
//   number.  I'm setting it to -1 in the arguments file (see
//   deconvolution_main.c for more detail). 
// 
// Subroutines: 
//      birth_death
//      mean_contribution
//      *mean_concentration
//      likelihood
//      *calc_death_rate_os
//      *calc_death_rate_strauss
//      calc_sr_strauss
// 
// Global variable definitions:
//      fitstart - The first time in hours that a pulse may occur
//      fitend   - The last time in hours that a pulse may occur
// 
//-----------------------------------------------------------------------------

// Include headers for files of used function definitions
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "birth_death.h"
#include "pulse_node.h" 
#include "calculations.h" 

// float.h included for EPS variable -- double check to see if this is
// working/necessary
//#include <float.h>
//#define EPS 1.0e-12

// Global variable definitions 
extern double fitstart;
extern double fitend;




//-----------------------------------------------------------------------------
//
// birth_death_strauss() 
//
//   this runs the birth-death step of the spatial BDMCMC process, using the
//   Strauss prior on pulse location.
//
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that exist
//     double **ts         - this is the matrix of observed data (a column
//                           of times and a column of log(concentration) 
//     Common_parms *parms - the current values of the common parameters
//     int N               - the number of observations in **ts
//     double *likeli      - the current value of the likelihood
//     int iter            - which iteration are we on
//     Priors *priors      - For priors->gamma, the repulsion parameter for
//                           Strauss process (hard core: gamma=0); and 
//                           priors->range, the range of repulsion (R) for
//                           Strauss process 
//
//   RETURNS:   
//     None                - all updates are made internally
//
//-----------------------------------------------------------------------------
void birth_death(Node_type *list, 
                 double **ts, 
                 Common_parms *parms, 
                 int N, 
                 double *likeli, 
                 //unsigned long *seed, 
                 int iter, 
                 int strauss,
                 Priors *priors) 
{

  // Declarations
  int i, j;                   // Generic counter
  int remove;                 // Ordered number of pulse to remove
  int num_node;               // Number of pulses
  int aaa = 0;                // Counter for # BD iterations
  int max_num_node = 60;      // Max number of pulses allowed before forced death
  double S = 0.0;             // Current progress in virtual time
  double Birth_rate;          // A constant birth rate needed for the process
  double T = 1.;              // Stop birthdeath loop when S exceeds T
  double full_likelihood;     // Current likelihood prior to birth or death
  double max;                 // Part of the calculation for Death_rate
  double Death_rate;          // Total death rate
  double *death_rate;         // Vector of death rates for each pulse
  double position;            // Value of new pulse's location
  double *partial_likelihood; // Array of likelihoods with pulse i removed
  //double *tmp;                // Used when drawing new pulse's mass and width
  Node_type *node;            // Pointer to current pulse
  Node_type *new_node;        // New pulse object
  // Spatial BD and Strauss declarations
  int sum_s_r = 0;            // Sum(S(R)) for birth of new pulse
  double birth_rate;          // 'Instantaneous' birth rate
  double papas_cif;           // Papangelou's cond'l intensity fn for birthed
  double b_ratio;             // Ratio of papas_cif/birth_rate for accept/rej
  double pulse_intensity;     // Prior intensity on pulse count poisson
  double new_theta;           // New theta draw for looping till > 0
  double new_eta;             // New eta draw for looping till > 0
  double new_tsd;             // New standard devation scaled by new_eta for
                              // drawing the new theta from a trunc t-distr.


  //-----------------------------------------
  // Prepare for birth death loop. Save and 
  // calculate values and allocate memory
  //-----------------------------------------
  full_likelihood = *likeli;
  Birth_rate = parms->nprior; 
  //tmp = (double *)calloc(2, sizeof(double));

  // If Strauss, calculate instantaneous birth rate and prior intensity
  if (strauss == 1) {
    birth_rate      = Birth_rate/(fitend-fitstart);
    pulse_intensity = parms->nprior/(fitend-fitstart);
  } 


  //-----------------------------------------
  // Start Birth-death loop 
  //   Run until break reached
  //-----------------------------------------
  while (1) {

    aaa++; // iters counter

    //-------------------------------------------
    // Calculate death rate for each component 
    // conditional on Strauss or order statistic
    //-------------------------------------------
    // Count number of pulses 
    num_node = 0;
    node = list->succ;
    while (node != NULL) {
      num_node++;
      node = node->succ;
    }

    // Calculate the likelihood if pulse i is removed
    i = 0;
    partial_likelihood = (double *)calloc(num_node, sizeof(double));
    node = list->succ;
    while (node != NULL) {
      partial_likelihood[i] = likelihood(list, ts, parms, N, node);
      i++;
      node = node->succ;
    }

    // Calculate death rate
    death_rate = NULL;
    if (strauss == 1) {
      death_rate = calc_death_rate_strauss(list, num_node, partial_likelihood,
                                           full_likelihood, birth_rate,
                                           pulse_intensity, priors);
    } else {
      death_rate = calc_death_rate_os(list, num_node, partial_likelihood,
                                      full_likelihood, Birth_rate, 
                                      parms->nprior, priors->orderstat);
    }


    //-------------------------------------------
    // Compute total death rate D
    //   Multi-step due to precision issues
    //-------------------------------------------
    if (death_rate != NULL) {

      Death_rate = death_rate[0];

      for (i = 1; i < num_node; i++) {
        max = (Death_rate > death_rate[i]) ? Death_rate : death_rate[i];
        Death_rate = log(exp(Death_rate - max) + 
                         exp(death_rate[i] - max)) + max;
      }

      for (i = 0; i < num_node; i++) {
        death_rate[i] -= Death_rate;
        death_rate[i]  = exp(death_rate[i]);
      }

      for (i = 1; i < num_node; i++) {
        death_rate[i] += death_rate[i-1];
      }

      if (Death_rate > 500 ) {
        Death_rate = 1e300;
      } else {
        Death_rate = exp(Death_rate); 
      }

    } else { 
      Death_rate = 0; 
    }

    if (num_node <= 1) { 
      Death_rate = 0; 
    }

    free(partial_likelihood);


    //-------------------------------------------
    // Calculate probability of birth
    //-------------------------------------------
    // Set Pr(Birth), force death/birth if too many/few
    if (num_node <= 1) { 
      max = 1.1;
    } else if (num_node >= max_num_node) { 
      max = -0.1;
    } else { 
      max = Birth_rate / (Birth_rate + Death_rate); 
    } 


    //-------------------------------------------
    // Update virtual time (how long to run BD step) 
    //   Draw from exp(B+D) and add to current S 
    //-------------------------------------------
    S += Rf_rexp(1/(Birth_rate + Death_rate)); //, seed);
    //Rprintf("S = %f; T = %f\n", S, T);
    // If S exceeds T or if we've run this too many times, break
    if (S > T)  { 
      //Rprintf("BD ran for %d iterations.\n", aaa);
      break;
    }
    if (aaa > 5000) {
      //Rprintf("BD ran for %d iterations.\n", aaa);
      break;
    } 


    //-------------------------------------------
    // Select Birth or Death and proceed with either
    //-------------------------------------------
    if (Rf_runif(0, 1) < max) { // If U < B/(B+D), a birth occurs 

      // Generate new position
      position = Rf_runif(fitstart, fitend);
      int accept_pos = 1;

      // If using Strauss prior, run accept/reject for new position
      if (strauss == 1) {
        sum_s_r    = calc_sr_strauss(position, list, list, priors);
        papas_cif  = pulse_intensity * pow(priors->gamma, sum_s_r);
        b_ratio    = papas_cif / birth_rate;

        accept_pos = (Rf_runif(0, 1) < b_ratio) ? 1 : 0;
      }

      // If it's a valid position, generate initial parms and insert node
      if (accept_pos == 1) {

        node = list->succ;
        // Initialize a new node
        new_node               = initialize_node();
        new_node->time         = position;

        for (j = 0; j < 2; j++) {
          new_theta = -1;
          new_eta = Rf_rgamma(2, 2);
          new_tsd = sqrt((parms->re_sd[j] * parms->re_sd[j]) / new_eta);
          while (new_theta < 0) {
            new_theta = Rf_rnorm(parms->theta[j], new_tsd);
          }
          new_node->theta[j] = new_theta;
          new_node->eta[j] = new_eta;
        }

        new_node->mean_contrib = (double *)calloc(N, sizeof(double));
        // Add it to the linklist and update values
        mean_contribution(new_node, ts, parms, N);
        insert_node(new_node, list);
        full_likelihood = likelihood(list, ts, parms, N, list);

      }

    } else { // Otherwise, a death occurs 

      // Pick a node to remove, find and remove it and update likelihood
      remove = one_rmultinom(death_rate, num_node) + 1; 
      node   = list;

      for (i = 0; i < remove; i++) { 
        node = node->succ; 
      }

      delete_node(node, list);
      full_likelihood = likelihood(list, ts, parms, N, list);

    } 

    // Clean up
    free(death_rate);

  } // End of BD while Loop 


  // Update likelihood before leaving Birthdeath.c
  *likeli = full_likelihood;

  free(death_rate);
  //free(tmp);

}




//-----------------------------------------------------------------------------
//
// mean_contribution() 
//
//   this updates a pulse's mean_contrib vector based on current values of
//   parameters
// 
//   ARGUMENTS: 
//     Node_type *node     - what pulse are we updating
//     double **ts         - this is the matrix of observed data (a column
//                           of times and a column of log(concentration)
//     Common_parms *parms - the current values of the common parameters
//     int N               - the number of observations in **ts
// 
//   RETURNS: 
//     None                - all updates are made internally
// 
//-----------------------------------------------------------------------------
void mean_contribution(Node_type *node, 
                       double **ts, 
                       Common_parms *parms, 
                       int N) {

  // Declarations --------------------------------
  int i;    // generic counter
  double x; // vars used in arithmetic
  double y; // vars used in arithmetic
  double z; // vars used in arithmetic
  double w; // vars used in arithmetic

  // Calculate mean contribution -----------------
  z  = node->theta[1] * parms->decay;
  y  = parms->decay * (0.5 * z  + node->time);
  z += node->time;
  w  = sqrt(2. * node->theta[1]);

  for (i = 0; i < N; i++) {
    x = (ts[i][0] - z) / w;
    x = Rf_pnorm5(x * sqrt(2), 0.0, 1.0, 1, 0);

    if (x == 0) {
      node->mean_contrib[i] = 0; 
    } else {
      node->mean_contrib[i] = node->theta[0] * x * 
                              exp(y - ts[i][0] * parms->decay);
    }
  }

  // Exit, returning nothing ---------------------

}




//-----------------------------------------------------------------------------
// 
// mean_concentration() 
//   this takes each pulse's mean_contrib vector and sums across them to get a
//   overall mean concentration at each time
// 
//   Arguments: 
//     Node_type *list     - this is the current list of pulses that exist
//     Common_parms *parms - the current values of the common parameters
//     int N               - the number of observations in **ts
//     Node_type *node_out - if we want, we can ignore a pulse
//     double **ts         - this is the matrix of observed data (a column
//                           of times and a column of log(concentration)
// 
//   Returns: 
//     x                   - the vector of sums
// 
//-----------------------------------------------------------------------------
double *mean_concentration(Node_type *list, Common_parms *parms, int N,
                           Node_type *node_out, double **ts) {

  // Declarations --------------------------------
  int i;           // Generic counter
  double *x;       // Vector of sums
  Node_type *node; // Pointer to current pulse

  // Allocated memory for concentration array ----
  x = (double *)calloc(N, sizeof(double));

  // Sum mean contribs at each time pt -----------
  node = list->succ;
  while (node != NULL) {
    if (node != node_out) {
      for (i = 0; i < N; i++) {
        x[i] += node->mean_contrib[i];
      }
    }
    node = node->succ;
  }

  // Add the baseline contribution ---------------
  for (i = 0; i < N; i++) {
    x[i] += parms->md[0];
    x[i] = log(x[i]);
  }

  // Exit returning concentration array ----------
  return x;

}




//-----------------------------------------------------------------------------
// 
// likelihood() 
//   computes the current likelihood using the observed log-concentrations
//   and mean concentration
// 
//   ARGUMENTS: 
//     Node_type *list     - this is the current list of pulses that exist
//     double **ts         - this is the matrix of observed data (a column
//                           of times and a column of log(concentration)
//     Common_parms *parms - the current values of the common parameters
//     int N               - the number of observations in **ts
//     Node_type *node_out - if we want, we can ignore a pulse
// 
//   RETURNS: 
//     x                   - the likelihood computed based on the inputs
//                           (scalar valued)
//  
//-----------------------------------------------------------------------------
double likelihood(Node_type *list, double **ts, Common_parms *parms, int N,
                  Node_type *node_out) {

  // Declarations --------------------------------
  int i;        // Generic counter
  double x = 0; // likelihood
  double *mean; // mean concentration

  // Sum across mean_contribs --------------------
  mean = mean_concentration(list, parms, N, node_out, ts);

  for (i = 0; i < N; i++) {
    x += (ts[i][1] - mean[i]) * (ts[i][1] - mean[i]);
  }

  x /= (-2.0 * parms->sigma);
  x += -0.5 * N * (1.8378771 + parms->lsigma);
  free(mean);

  // Exit, returning likelihood scalar -----------
  return x;

}




//-----------------------------------------------------------------------------
// 
// calc_death_rate_os() 
//   calculates a vector of death rates, one for each existing pulse
//   
//   ARGUMENTS: 
//     Node_type *list            - Current linklist of pulses 
//     int num_node               - Current number of pulses
//     double *partial_likelihood - Vector of partial likelihoods, where the
//                                  ith element represents the likelihood with
//                                  the ith pulse removed
//     double full_likelihood     - Value of the full likelihood
//     double Birth_rate          - Value of the birth rate
//     double r                   - Prior on pulse count (parms->nprior)
//     int mmm                    - Order stat to use for prior on location
//                                  (priors->orderstat)
// 
//   RETURNS:
//     death_rate                 - Vector where the ith element represents 
//                                  the death rate of the ith pulse
// 
//-----------------------------------------------------------------------------
double *calc_death_rate_os(Node_type *list, 
                           int num_node, 
                           double *partial_likelihood, 
                           double full_likelihood, 
                           double Birth_rate, 
                           double r,
                           int mmm) {

  int i;              // Generic counter
  int j;              // Generic counter
  double x;           // Individual death rate
  double *death_rate; // Vector of death rates
  double coef_denom;  // Part of death rate calculation
  double coef_num;    // Part of death rate calculation
  Node_type *node;    // Pointer to node in linklist

  if (num_node > 1) {
    death_rate = (double *)calloc(num_node, sizeof(double));
    node = list->succ;
    i = 0;

    // Calculate the coefficient of the distribution of the taus conditional
    // on number of pulses. In this portion, I have an extra num_node in the
    // numerator 
    coef_num = 1;
    for (j = 1; j < mmm; j++) { coef_num *= j; }
    coef_denom = mmm * num_node;
    for (j = 1; j < mmm; j++) { coef_denom *= (mmm * num_node + j); }

    while (node != NULL) {

      if (mmm > 0) {
        // This computes the portion of death rate due to the Poisson prior on
        // N, the birth rate, and the likelihood */
        // num_node = number of pulses
        // r = prior on pulse count
        x = log(num_node * Birth_rate / r) + 
            partial_likelihood[i] - full_likelihood;
      } else {
        x = log(Birth_rate / num_node) + 
            partial_likelihood[i] - full_likelihood;
      }

      // Now, compute the portion of the death rate due to distribution of the
      // taus conditional on number of pulses.
      if (mmm > 1) {

        if (i ==  0) {

          // If we are on the first pulse
          x += (mmm-1) * log(fitend - fitstart) +
               (mmm-1) * log((node->succ->time - fitstart) /
                             ((node->succ->time - node->time) * 
                              (node->time - fitstart))) +
               log(coef_num / coef_denom);

        } else if (i > 0) {

          // If we are not on the first pulse
          if (node->succ) {
            x += (mmm-1) * log(fitend - fitstart) +
                 (mmm-1) * log((node->succ->time - node->pred->time) /
                               ((node->succ->time - node->time) * 
                                (node->time - node->pred->time))) +
                 log(coef_num / coef_denom);

          } else { // If we are not on the first or last pulse 
            x += (mmm-1) * log(fitend - fitstart) + 
                 (mmm-1) * log((fitend - node->pred->time) / 
                               ((fitend - node->time) * 
                                (node->time - node->pred->time))) +
                 log(coef_num / coef_denom);
          }

        }

      }

      // if pulse equal to NaN, then set to large value/guarantee death
      // necessary when starting value is < fitstart (orderstat part of
      // death calc causes this due to neg. value in log) 
      if (isnan(x)) {
        x = 1e300;
      }
      // Save to death rate vector 
      death_rate[i] = x;

      // Advance to next pulse
      node = node->succ;

      // Advance counter 1
      i++;

    }

  } else {

    // If we have 0/1 pulses return NULL/0
    if (num_node == 0) {
      return NULL;
    } else {
      death_rate = (double *)calloc(num_node,sizeof(double));
      death_rate[0] = -1e300;
    }

  }

  return death_rate;  

}




//-----------------------------------------------------------------------------
// 
//  calc_death_rate_strauss()
// 
//    Calculates a vector of death rates, one for each existing pulse
//    
//    Arguments: 
//      Node_type *list            - this is the current list of pulses that
//                                   exist
//      int num_node               - current number of pulses
//      double *partial_likelihood - vector of likelihoods, where the ith
//                                   element represents the likelihood with
//                                   the ith pulse removed
//      double full_likelihood     - value of the full likelihood
//      double Birth_rate          - value of the birth rate
// 
//    Returns:
//      death_rate                 - a vector where the ith element
//                                   represents the death rate of the ith
//                                   pulse
// Notes: 
//   Two possible version:
//      1) if we can use distr of birth_rate, it's simple -- just partial -
//      full likelihoods.
//      2) if we have to use our set birth_rate (which I suspect based on
//      Stephens2000), we have to multiply by papas^-1 * birth_rate
//      
//   Resolution:
//     1) is the correct approach. We aren't actually setting the birth rate.
//     Instead we are providing the intensity parameter as the strauss and as
//     long as death rate and birth rate intensity parms are the same, they
//     cancel. Option 2) would require not using an accept/reject in birth step
//     and instead generating from the actual Strauss density (i.e. really
//     difficult and computationally expensive.
// 
//-----------------------------------------------------------------------------
double *calc_death_rate_strauss(Node_type *list, 
                                int num_node, 
                                double *partial_likelihood, 
                                double full_likelihood, 
                                double birth_rate,
                                double pulse_intensity, 
                                Priors *priors) {

  // Declarations --------------------------------
  int i = 0;          // Generic counter
  double x;           // Variable for death rate
  double *death_rate; // Individual death rate array;
  Node_type *node;    // Pointer to pulses linked list

  // Begin death rate calculations ---------------
  if (num_node > 1) {

    death_rate = (double *)calloc(num_node, sizeof(double));
    node = list->succ;

    // Calculate death rates (priors all cancel), save, 
    // and increment pulse and counter
    while (node != NULL) {
      x = partial_likelihood[i] - full_likelihood; 
      death_rate[i] = x;
      node = node->succ;
      i++;
    }

  } else { 

    // if we have 0 or 1 pulses return NULL/0 
    // (i.e. don't kill anything)
    if (num_node==0) { 
      return NULL;
    } else {
      death_rate = (double *)calloc(num_node, sizeof(double));
      death_rate[0] = -1e300;
    }

  }

  // Exit, returning death rate array ------------
  return death_rate;  

}




//-----------------------------------------------------------------------------
// 
//  calc_sr_strauss()
// 
//    Calculates sum(S(R)), the exponent on the gamma parameter in the Strauss
//    process/prior for pulse location. Used for Strauss prior 
//    (prior->gamma >= 0) prior in birth_death and mh_time_strauss.
// 
//    Arguments: 
//      position    position to test nodes against (current or proposal,
//                  depending on context)
//      *list       linked list of all pulse nodes
//      *node_out   optional node to exclude
//      *priors     prior parameters
// 
//    Returns: 
//      sum(S(R))   scalar value for sum of # pulses too close to each other
// 
//-----------------------------------------------------------------------------
int calc_sr_strauss(double position, Node_type *list, Node_type *node_out,
                    Priors *priors) {

  int s_r;           // Sum of indicators where diff < 20
  double difference; // Time difference
  Node_type *node;   // Pointer to node for comparing to 'position'

  s_r = 0;
  node = list->succ;

  while(node != NULL) {
    if (node != node_out) {
      // skip if node is same that position is from;
      difference = fabs(position - node->time);
      // increment by 1 if diff<R
      s_r = (difference < priors->range) ? s_r + 1 : s_r; 
    }
    node = node->succ;
  }

  return(s_r); 

}




//------------------------------------------------------------------------------
// END OF FILE
//------------------------------------------------------------------------------

