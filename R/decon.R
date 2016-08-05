#-------------------------------------------------------------------------------
# decon.R - Wrapper for fitting decon
#-------------------------------------------------------------------------------

#decon <- function(data, priors, starting.vals, proposal.vars, seeds) {
#
#  .Call(decon_, data, priors, starting.vals, proposal.vars, seeds)
#
#}


#' test_inout
#'
#' testing reading in data, converting to c data types and converting back.
#'   test list/dataframe to 
#' 
#' @param data named vector of parameters
#' @keywords pulse simulation
#' @export
#' @useDynLib pulsatile testc
#' @useDynLib pulsatile testspec
#' @examples none currently
#' test_inout()
test_inout <- function(x, ...) {
  
  if ("pulse_spec" %in% class(x)) { 
    .Call(testspec, x, PACKAGE = "pulsatile")

  } else if (is.data.frame(x)) {
    .Call(testc, x, PACKAGE = "pulsatile")
  }

}


#' show_args
#'
#' Example via R source code: http://git.io/v0zFx
#' 
#' @useDynLib pulsatile showArgs1
#' @param data named vector of parameters
#' @keywords pulse simulation
#' @export
#' @examples none currently
#' show_args()
show_args <- function(.data) {
  
  .Call(showArgs1, .data, PACKAGE = "pulsatile")

}


#' fit_pulse
#' 
#' Primary function for fitting deconvulation model for time series of pulsatile
#' hormone data
#' 
#' @useDynLib pulsatile decon
#' @param pulse_spec_obj 
#' @keywords pulse fit
#' @export
#' @examples none currently
#' fit_pulse()
fit_pulse <- function(pulse_spec_obj) {

  .Call(decon, PACKAGE = "pulsatile",
        pulse_spec_obj$model$data,
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
        pulse_spec_obj$prior$error$beta,
        pulse_spec_obj$starting_values$pulse_mass$mean,
        pulse_spec_obj$starting_values$pulse_mass$sd,
        pulse_spec_obj$starting_values$pulse_width$mean,
        pulse_spec_obj$starting_values$pulse_width$sd,
        pulse_spec_obj$starting_values$baseline$mean,
        pulse_spec_obj$starting_values$halflife$mean,
        pulse_spec_obj$starting_values$error$var,
        pulse_spec_obj$proposal_variances$mean_pulse_mass,
        pulse_spec_obj$proposal_variances$mean_pulse_width,
        pulse_spec_obj$proposal_variances$pulse_mass,
        pulse_spec_obj$proposal_variances$pulse_width,
        pulse_spec_obj$proposal_variances$pulse_location,
        pulse_spec_obj$proposal_variances$baseline,
        pulse_spec_obj$proposal_variances$halflife)

}
