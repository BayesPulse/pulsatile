///-----------------------------------------------------------------------------
///
/// FILE: decon.h
///
/// DESCRIPTION: 
///   Function definitions for decon.c
/// 
///-----------------------------------------------------------------------------

#include <R.h>
#include <Rinternals.h>

// Function Definitions
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
           SEXP pv_halflife);
double **convert_data(SEXP indata);
SEXP getListElement(SEXP list, const char *str);
