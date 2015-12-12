///-----------------------------------------------------------------------------
///
/// FILE: birthdeath_strauss.h 
///
/// DESCRIPTION: 
///   Function definitions for birthdeath_strauss.c
/// 
///-----------------------------------------------------------------------------

#ifndef BIRTHDEATH_STRAUSS_H
#define BIRTHDEATH_STRAUSS_H

// Include this here since defn of structures used in arguments
#include "deconvolution_main.h"


// Birth-death algorithm for poisson/order-statistic or Strauss priors on pulse
// count/location
void birth_death(Node_type *list,  //
                 double **ts, 
                 Common_parms *parms, 
                 int N, 
                 double *likeli, 
                 unsigned long *seed, 
                 int iter, 
                 Priors *priors);

// Calculates mean contribution for a given pulse at each time point
void mean_contribution(Node_type *node, 
                       double **ts, 
                       Common_parms *parms, 
                       int N);

// Calculates mean concentration by summing all pulse's mean contribution at
// each time point
double *mean_concentration(Node_type *list, 
                           Common_parms *parms, 
                           int N,
                           Node_type *node_out, 
                           double **ts);

// Calculates likelihood for current set of pulses and parameter estimates
double likelihood(Node_type *list, 
                  double **ts,
                  Common_parms *parms, 
                  int N,
                  Node_type *node_out);

// Calculates death rate of each pulse in provided linklist for the order
// statistic prior on pulse location
double *calc_death_rate_os(Node_type *list, 
                           int num_node, 
                           double *partial_likelihood, 
                           double full_likelihood, 
                           double Birth_rate, 
                           double r);

// Calculates death rate of each pulse in provided linklist for the Strauss
// prior on pulse count and location
double *calc_death_rate_strauss(Node_type *list, 
                                int num_node, 
                                double *partial_likelihood, 
                                double full_likelihood, 
                                double birth_rate,
                                double pulse_intensity, 
                                Priors *priors);

// Calculates sum[S(R)] portion of the Strauss prior.  sum[S(R)] is the count of
// pulses that are less than R minutes from each other. 
int calc_sr_strauss(double position, 
                    Node_type *list, 
                    Node_type *node_out,
                    Priors *priors);

#endif // BIRTHDEATH_STRAUSS_H

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

