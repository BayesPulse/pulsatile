
#-------------------------------------------------------------------------------
# Functions for creating a pulse model specification
#-------------------------------------------------------------------------------

#' pulse_spec.R
#'
#' Generates a pulse_spec object -- the specification object required for
#' fitting a fit_pulse model.
#'   
#' 
#' @param location_prior_type Takes on two values: "order-statistic" and
#' "strauss". "order-statistic" uses every third order statistic of a Uniform
#' distribution for the pulse location prior and requires specification of the
#' prior parameter for mean pulse count ("prior_mean_pulse_count").
#' "strauss" uses the Strauss interacting point-process as a prior and requires
#' specification of "prior_mean_pulse_count", "prior_location_gamma", and
#' "prior_location_range".
#' @param prior_mass_mean mass mean hyperparm
#' @param prior_mass_var mass variance hyperparm
#' @param prior_width_mean width mean hyperparm
#' @param prior_width_var width variance hyperparm
#' @param prior_baseline_mean mean of prior on baseline
#' @param prior_baseline_var variance of prior on baseline
#' @param prior_halflife_mean mean of prior on half-life
#' @param prior_halflife_var variance of prior on half-life
#' @param prior_error_alpha placeholder
#' @param prior_error_beta placeholder
#' @param prior_location_gamma placeholder
#' @param prior_location_range placeholder
#' @param prior_max_sd_mass placeholder
#' @param prior_max_sd_width placeholder
#' @param prior_mean_pulse_count placeholder
#' @param sv_mass_mean placeholder
#' @param sv_width_mean placeholder
#' @param sv_baseline_mean placeholder
#' @param sv_halflife_mean placeholder
#' @param sv_error_var placeholder
#' @param sv_mass_sd placeholder
#' @param sv_width_sd placeholder
#' @param pv_baseline placeholder
#' @param pv_halflife placeholder
#' @param pv_mean_pulse_mass placeholder
#' @param pv_mean_pulse_width placeholder
#' @param pv_pulse_mass placeholder
#' @param pv_pulse_width placeholder
#' @param pv_pulse_location placeholder
#' @export
#' @keywords pulse simulation
#' pulse_spec()
pulse_spec <-
  function(location_prior_type = "order-statistic",
           prior_mass_mean        = 1.50,
           prior_mass_var         = 10,
           prior_width_mean       = 3.5,
           prior_width_var        = 10,
           prior_baseline_mean    = 2.6,
           prior_baseline_var     = 100,
           prior_halflife_mean    = 45,
           prior_halflife_var     = 100,
           prior_error_alpha      = 0.0001,
           prior_error_beta       = 0.0001,
           prior_location_gamma   = NULL,
           prior_location_range   = NULL,
           prior_max_sd_mass      = 10,
           prior_max_sd_width     = 10,
           prior_mean_pulse_count = 12,
           sv_mass_mean           = 1.5,
           sv_width_mean          = 3.5,
           sv_baseline_mean       = 2.6,
           sv_halflife_mean       = 45,
           sv_error_var           = 0.25,
           sv_mass_sd             = 1,
           sv_width_sd            = 1,
           pv_baseline            = 0.5,
           pv_halflife            = 45,
           pv_mean_pulse_mass     = 2,
           pv_mean_pulse_width    = 2,
           pv_pulse_mass          = 2,
           pv_pulse_width         = 50,
           pv_pulse_location      = 10) 
  {

    # Structure for single-subject, strauss prior model.
    # NOTE: sv's use std dev, while priors use variances (want this consistent
    # for any reason?)
    # NOTE: need more clear label for max_sd's 
    ps_obj <- 
      structure(
        list(location_prior_type = location_prior_type,
             priors = list(pulse_mass     = list(mean  = prior_mass_mean,
                                                 var   = prior_mass_var),
                           pulse_width    = list(mean  = prior_width_mean,
                                                 var   = prior_width_var),
                           pulse_location = list(gamma = prior_location_gamma,
                                                 range = prior_location_range,
                                                 count = prior_mean_pulse_count),
                           max_sd         = list(mass  = prior_max_sd_mass,
                                                 width = prior_max_sd_width),
                           baseline       = list(mean  = prior_baseline_mean,
                                                 var   = prior_baseline_var),
                           halflife       = list(mean  = prior_halflife_mean,
                                                 var   = prior_halflife_var),
                           error          = list(alpha = prior_error_alpha,
                                                 beta  = prior_error_beta)),
             starting_values = list(pulse_mass  = list(mean = sv_mass_mean,
                                                       sd   = sv_mass_sd),
                                    pulse_width = list(mean = sv_width_mean,
                                                       sd   = sv_width_sd),
                                    baseline    = list(mean = sv_baseline_mean),
                                    halflife    = list(mean = sv_halflife_mean),
                                    error       = list(var  = sv_error_var)),
             proposal_variances = list(mean_pulse_mass  = pv_mean_pulse_mass,
                                       mean_pulse_width = pv_mean_pulse_width,
                                       pulse_mass       = pv_pulse_mass,
                                       pulse_width      = pv_pulse_width,
                                       pulse_location   = pv_pulse_location,
                                       baseline         = pv_baseline,
                                       halflife         = pv_halflife)),
        class = "pulse_spec")

    return(ps_obj)

  }


#' @export 
print.pulse_spec <- function(x, ...) {

  cat("\nBayesian time-series analysis of pulsatile hormone data: 
      Model Specification Object\n\n")
  cat("Model type:", paste0(x$model$model, "\n"))
  cat("Number of iterations:", 
      formatC(x$model$iterations, format = "d", big.mark = ","), "\n")
  cat("\nPulse mass:\n")
  cat("   prior mean =", x$priors$pulse_mass$mean, "\n") 
  cat("   prior variance =", x$priors$pulse_mass$var, "\n") 
  cat("   mean starting value =", x$starting_values$pulse_mass$mean, "\n") 
  cat("   SD starting value =", x$starting_values$pulse_mass$sd, "\n") 
  cat("   proposal variance =", x$proposal_variances$pulse_mass, "\n")
  cat("   mean proposal variance =", x$proposal_variances$mean_pulse_mass, "\n")

}


################################################################################
# End of file # End of file  # End of file  # End of file  # End of file  # End
################################################################################
