///-----------------------------------------------------------------------------
///
/// FILE: mcmc.h
///
/// DESCRIPTION: 
///   Function definitions for linklistv3.c
/// 
///-----------------------------------------------------------------------------

#ifndef MCMC_H
#define MCMC_H

// Include this here since defn of structures used in arguments

#endif 

#include "r_interface.h"
#include "pulse_node.h"

void mcmc(Node_type *list, Common_parms *parms, double **ts, long iter, int N,
          int NN, int strauss, Priors *priors, SEXP common, SEXP pulse_chains,
          double propsd[]);

void mh_time_strauss(Node_type *list, Common_parms *parms, double **ts, 
                     double *like, int N, double sd, Priors *priors, 
                     long *atime, long *ntime);

void mh_time_os(Node_type *list, Common_parms *parms, double **ts, 
                double *like, int N, double v, long *atime, long *ntime);

void mh_mu_delta(Node_type *list, Common_parms *parms, Priors *priors, 
                 double **ts, double *like, int N, int num_node, double **var,
                 long *adelta, long *ndelta);

void draw_fixed_effects(Node_type *list, Priors *priors, Common_parms *parms);

void draw_re_sd(Node_type *list, Priors *priors, Common_parms *parms, 
                double v1, double v2, long *arevm, long *nrevm, 
                long *arevw, long *nrevw);

void draw_random_effects(double **ts, Node_type *list, Common_parms *parms, 
                         int N, double *like, double v1,  double v2, 
                         long *arem, long *nrem, long *arew, 
                         long *nrew);

double error_squared(double **ts, Node_type *list, Common_parms *parms, int N);

void adjust_acceptance(double x, double *X);

void adjust2_acceptance(double x, double **X, double corr);

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

