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

#include "decon_test.h"
#include "pulse_node.h"

void mcmc(Node_type *list, Common_parms *parms, double **ts, long iter, int N,
          int NN, Priors *priors, SEXP common, SEXP parm, double propsd[]);

void mh_time_strauss(Node_type *list, Common_parms *parms, double **ts, 
                     double *like, int N, double v, Priors *priors);

void mh_time_os(Node_type *list, Common_parms *parms, double **ts, 
                double *like, int N, double v);

void mh_mu_delta(Node_type *list, Common_parms *parms, Priors *priors, 
                 double **ts, double *like, int N, int num_node, double **var);

void draw_fixed_effects(Node_type *list, Priors *priors, Common_parms *parms);

void draw_re_sd(Node_type *list, Priors *priors, Common_parms *parms, 
                double v1, double v2 );

void draw_random_effects(double **ts, Node_type *list, Common_parms *parms, 
                         int N, double *like, double v1,  double v2 );

double error_squared(double **ts, Node_type *list, Common_parms *parms, int N);

void adjust_acceptance(double x, double *X);

void adjust2_acceptance(double x, double **X, double corr);

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

