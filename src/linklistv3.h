///-----------------------------------------------------------------------------
///
/// FILE: linklistv3.h
///
/// DESCRIPTION: 
///   Function definitions for linklistv3.c
/// 
///-----------------------------------------------------------------------------

#ifndef LINKLISTV3_H
#define LINKLISTV3_H

// Include this here since defn of structures used in arguments
#include "deconvolution_main.h"

void mcmc(Node_type *list, 
          Common_parms *parms, 
          double **ts, 
          long iter, 
          int N,
          Priors *priors, 
          unsigned long *seed, 
          char *file1, 
          char *file2, 
          double propvar[]);

void mh_time_strauss(Node_type *list, 
                     Common_parms *parms, 
                     double **ts, 
                     double *like, 
                     int N, 
                     unsigned long *seed, 
                     double v, 
                     Priors *priors);

void mh_time_os(Node_type *list, 
                Common_parms *parms, 
                double **ts, 
                double *like,
                int N, 
                unsigned long *seed, 
                double v);

void mh_mu_delta(Node_type *list, 
                 Common_parms *parms, 
                 Priors *priors, 
                 double **ts, 
                 double *like, 
                 int N, 
                 int num_node, 
                 unsigned long *seed, 
                 double **var);

void draw_fixed_effects(Node_type *list, 
                        Priors *priors, 
                        Common_parms *parms,
                        unsigned long *seed);

void draw_re_sd(Node_type *list, 
                Priors *priors, 
                Common_parms *parms,
                double v1, 
                double v2, 
                unsigned long *seed);

void draw_random_effects(double **ts, 
                         Node_type *list, 
                         Common_parms *parms, 
                         int N, 
                         double *like, 
                         double v1,  
                         double v2, 
                         unsigned long *seed);

double error_squared(double **ts, 
                     Node_type *list, 
                     Common_parms *parms, 
                     int N);

void adjust_acceptance(double x, 
                       double *X);

void adjust2_acceptance(double x, 
                        double **X, 
                        double corr);






#endif // LINKLISTV3_H

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

