#-------------------------------------------------------------------------------
# decon.R - Wrapper for fitting decon
#-------------------------------------------------------------------------------


#' fit_pulse
#' 
#' Primary function for fitting deconvulation model for time series of pulsatile
#' hormone data
#' 
#' @useDynLib pulsatile decon_input
#' @param pulse_spec_obj 
#' @keywords pulse fit
#' @export
#' @examples none currently
#' fit_pulse()
fit_pulse <- function(pulse_spec_obj) {

  .Call(decon_input, PACKAGE = "pulsatile",
        as.data.frame(pulse_spec_obj$model$data),
        pulse_spec_obj$model$model,
        as.integer(pulse_spec_obj$model$thin),
        as.integer(pulse_spec_obj$model$iterations),
        pulse_spec_obj$priors$pulse_mass$mean,
        pulse_spec_obj$priors$pulse_mass$var,
        pulse_spec_obj$priors$pulse_width$mean,
        pulse_spec_obj$priors$pulse_width$var,
        pulse_spec_obj$priors$pulse_location$gamma,
        pulse_spec_obj$priors$pulse_location$range,
        pulse_spec_obj$priors$pulse_location$count,
        pulse_spec_obj$priors$max_sd$mass,
        pulse_spec_obj$priors$max_sd$width,
        pulse_spec_obj$prior$baseline$mean,
        pulse_spec_obj$prior$baseline$var,
        pulse_spec_obj$prior$halflife$mean,
        pulse_spec_obj$prior$halflife$var,
        pulse_spec_obj$prior$error$alpha,
        pulse_spec_obj$prior$error$beta) #,
        #pulse_spec_obj$starting_values$pulse_mass$mean,
        #pulse_spec_obj$starting_values$pulse_mass$sd,
        #pulse_spec_obj$starting_values$pulse_width$mean,
        #pulse_spec_obj$starting_values$pulse_width$sd,
        #pulse_spec_obj$starting_values$baseline$mean,
        #pulse_spec_obj$starting_values$halflife$mean,
        #pulse_spec_obj$starting_values$error$var,
        #pulse_spec_obj$proposal_variances$mean_pulse_mass,
        #pulse_spec_obj$proposal_variances$mean_pulse_width,
        #pulse_spec_obj$proposal_variances$pulse_mass,
        #pulse_spec_obj$proposal_variances$pulse_width,
        #pulse_spec_obj$proposal_variances$pulse_location,
        #pulse_spec_obj$proposal_variances$baseline,
        #pulse_spec_obj$proposal_variances$halflife)

}
