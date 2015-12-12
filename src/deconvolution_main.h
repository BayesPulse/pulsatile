///-----------------------------------------------------------------------------
///
/// FILE: deconvolution_main.h
///
/// DESCRIPTION: 
///   Data structure definitions for entire program.
/// 
///-----------------------------------------------------------------------------

#ifndef DECONVOLUTION_MAIN_H
#define DECONVOLUTION_MAIN_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

typedef struct node_tag {
    struct node_tag *succ; // next node
    struct node_tag *pred; // previous node
    double time;           // pulse location for individual pulses
    double theta[2];       // 0: individual pulse mass; 1: individual pulse variance
    double *mean_contrib;  // vector of secretion contribution for this pulse at
                           // each timepoint
} Node_type;


typedef struct {
    double md[2];    // 0: baseline; 1: half-life
    double sigma;    // model error (variance)
    double lsigma;   // log of model error (may not be used)
    double theta[2]; // mean pulse mass and width on log scale
    double *re_sd;   // variance/covariance matrix for theta
    double decay;    // decay rate converted from above half-life
    double nprior;   // prior number of pulses
} Common_parms;

typedef struct {
    double meanbh[2];    // prior mean on baseline and halflife
    double varbh[2];     // variance on prior of baseline and halflife
    double *fe_mean;     // array of prior mean for mean pulse mass and width
    double *fe_variance; // variance of prior for mean pulse mass
    double *re_sdmax;    // 'b' on uniform prior for random effects standard deviation
    double alpha;        // alpha: prior parameters for model error (inverse gamma)
    double beta;         // beta : ""
    double gamma;        // Gamma, repulsion parameter for Strauss prior
    double range;        // Range parameter for Strauss prior
} Priors;


#endif // DECONVOLUTION_MAIN_H

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

