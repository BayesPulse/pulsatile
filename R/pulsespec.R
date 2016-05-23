#-------------------------------------------------------------------------------
# Functions for creating a pulse model specification
#-------------------------------------------------------------------------------

#' pulse_spec.R
#'
#' Generates a pulse_spec object -- the specification object required for
#' fitting a fit_pulse model.
#'   
#' 
#' @param iterations Number of iterations to run MCMC_
#' @param model Model type to fit. One of "single-series", "population",
#' "single-series associational", "population associational".
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
#' @param max_sd_mass placeholder
#' @param max_sd_width placeholder
#' @param prior_mean_num_pulses placeholder
#' @param sv_mass_mean placeholder
#' @param sv_width_mean placeholder
#' @param sv_baseline_mean placeholder
#' @param sv_halflife_mean placeholder
#' @param sv_error_var placeholder
#' @param sv_mass_sd placeholder
#' @param sv_width_sd placeholder
#' @param pv_baseline placeholder
#' @param pv_halflife placeholder
#' @param pv_mass placeholder
#' @param pv_width placeholder
#' @param pv_re_mass placeholder
#' @param pv_re_width placeholder
#' @param pv_pulselocations placeholder
#' @keywords pulse simulation
#' @export
#' @examples none currently
#' @seealso \code{\link{summary.pulse_spec}} \code{\link{print.pulse_spec}} 
#' pulse_spec()
pulse_spec <-
  function(.data = NULL,
           iterations = 250000,
           model = c("single-series"), #, "population", "single-series
                     #associational", "population associational"),
           prior_mass_mean       = 1.50,
           prior_mass_var        = 10,
           prior_width_mean      = 3.5,
           prior_width_var       = 10,
           prior_baseline_mean   = 2.6,
           prior_baseline_var    = 100,
           prior_halflife_mean   = 45,
           prior_halflife_var    = 100,
           prior_error_alpha     = 0.0001,
           prior_error_beta      = 0.0001,
           prior_location_gamma  = 0,
           prior_location_range  = 40,
           prior_max_sd_mass     = 10,
           prior_max_sd_width    = 10,
           prior_mean_num_pulses = 12,
           sv_mass_mean          = 1.5,
           sv_width_mean         = 3.5,
           sv_baseline_mean      = 2.6,
           sv_halflife_mean      = 45,
           sv_error_var          = 0.25,
           sv_mass_sd            = 1,
           sv_width.sd           = 1,
           pv_baseline           = 0.5,
           pv_halflife           = 45,
           pv_mass               = 2,
           pv_width              = 2,
           pv_re_mass            = 2,
           pv_re_width           = 50,
           pv_pulselocations     = 10) 
  {

    ps_obj <- 
      structure(list(model = list("model" = model, iterations = iterations, data = .data),
                     priors = list(pulse_mass = 
                                   list(mean = prior_mass_mean,
                                        var  = prior_mass_var),
                                   pulse_width = 
                                     list(mean = prior_width_mean,
                                          var  = prior_width_var),
                                   baseline_concentration = 
                                     list(mean = prior_baseline_mean,
                                          var  = prior_baseline_var),
                                   halflife =
                                     list(mean = prior_halflife_mean,
                                          var  = prior_halflife_var)),
                     starting_values = 
                       list(),
                     proposal_variances =
                       list()),
                class = "pulse_spec")



    #class(ps_obj) <- c("pulse_spec")
    return(ps_obj)

  }


print <- function(x) { useMethod("pulse_spec") }

summary <- function(x) { useMethod("pulse_spec") }

print.pulse_spec <- function(x) {
  cat("\nModel type: ", paste0(x$model$model, "\n\n"))

}



################################################################################
# End of file # End of file  # End of file  # End of file  # End of file  # End
################################################################################
