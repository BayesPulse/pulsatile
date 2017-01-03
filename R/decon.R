#-------------------------------------------------------------------------------
# decon.R - Wrapper for fitting decon
#-------------------------------------------------------------------------------


#' fit_pulse
#' 
#' Primary function for fitting deconvulation model for time series of pulsatile
#' hormone data
#' 
#' @useDynLib pulsatile decon_r_interface
#' @param formula A formula of the form 'observed concentration ~ time'.
#' @param data A dataset containing the variables in formula.
#' @param iterations Number of iterations for MCMC
#' @param thin Thinning to apply to MCMC chains (i.e. Keep every 'thin'th
#' sample).
#' @param pulse_spec_obj An object of class 'pulse_spec_obj', created by
#' pulse_spec(), specifying the priors, starting values, and proposal variances
#' to use.
#' @param use_tibble Return chains as tbl_df class data frames, from the tibble
#' package.  Mostly used for the print.tbl_df method, which limits the rows and
#' columns printed to those which fit in the console.
#' @import tibble
#' @export
#' @keywords pulse fit
#' fit_pulse()
fit_pulse <- function(formula,
                      data,
                      iterations = 2500,
                      #model_type = c("single-series"), #, "population", "single-series
                      #associational", "population associational"),
                      thin       = 50,
                      pulse_spec_obj = pulse_spec(),
                      use_tibble = TRUE) {

  #model <- match.arg(model)
  # ideas via survival::coxph 
  Call  <- match.call()
  arg_indx <- match(c("formula", "data", "iterations", "thin", "pulse_spec_obj"), 
                    names(Call), nomatch = 0)

  if (class(pulse_spec_obj) != "pulse_spec") {
    stop("pulse_spec_obj is invalid -- see the fit_pulse() and pulse_spec()
         documentation.")
  }

  # 
  #this_env <- environment()
  #list2env(pulse_spec_obj, envir = this_env)

  rtn <- .Call(decon_r_interface, PACKAGE = "pulsatile",
               data,
               "single-series",
               as.integer(thin),
               as.integer(iterations),
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

  common_chain <- as.data.frame(rtn[[1]])
  pulse_chain  <- as.data.frame(do.call(rbind, rtn[[2]]))
  if (use_tibble) {
    common_chain = tibble::as_data_frame(common_chain)
    pulse_chain = tibble::as_data_frame(pulse_chain)
  }

  rtn_obj <- 
    structure(list("model"        = "single-subject",
                   "call"         = Call,
                   "common_chain" = common_chain,
                   "pulse_chain"  = pulse_chain,
                   "data"         = data,
                   "options"      = list("thinning"   = thin, 
                                         "iterations" = iterations),
                   "spec"         = pulse_spec_obj),
              class = "pulse_fit")

  return(rtn_obj)

}
