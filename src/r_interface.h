//
// FILE: r_interface.h
// DESCRIPTION: Definitions for r_interface.c
// 

#define R_NO_REMAP 
#include <R.h>
#include <Rinternals.h>

#ifndef DECON_TEST_H
#define DECON_TEST_H

typedef struct {
    double meanbh[2];    // prior mean on baseline and halflife
    double varbh[2];     // variance on prior of baseline and halflife
    double *fe_mean;     // array of prior mean for mean pulse mass and width
    double *fe_variance; // variance of prior for mean pulse mass
    double *re_sdmax;    // 'b' on uniform prior for random effects standard deviation
    double err_alpha;    // alpha: prior parameters for model error (inverse gamma)
    double err_beta;     // beta : ""
    double orderstat;    // # orderstat to use for uniform prior on location
    double gamma;        // Gamma, repulsion parameter for Strauss prior
    double range;        // Range parameter for Strauss prior
} Priors;

typedef struct {
    double md[2];    // 0: baseline; 1: half-life
    double sigma;    // model error (variance)
    double lsigma;   // log of model error (may not be used)
    double theta[2]; // mean pulse mass and width on log scale (TODO: converting this to natural scale for truncated t-distr)
    double *re_sd;   // variance/covariance matrix for theta
    double decay;    // decay rate converted from above half-life
    double nprior;   // prior number of pulses
} Common_parms;

SEXP decon_r_interface(SEXP indata,
                       SEXP model,
                       SEXP thin,
                       //SEXP burnin,
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
                       );
double **convert_data(SEXP indata, int nrow); 
void deallocate_data(double **data, int nrow);
SEXP getListElement(SEXP list, const char *str);
//#endif // DECONVOLUTION_MAIN_H

#endif

//------------------------------------------------------------------------------
// End of file
//------------------------------------------------------------------------------

