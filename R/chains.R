#   Extract heirarchical chains
#     - hormone1_chain() ??? driver_chain()
#     - hormone2_chain() ??? response_chain()
#     - pulse_chain()
#     - subject_chain()
#     - pop_chain()
#     - chains() --> list/str chains available and functions to extract them
#

# #' fit_pulse
# #' 
# #' Primary function for fitting deconvolution model for time series of pulsatile
# #' hormone data
# #' 
# #' @useDynLib pulsatile decon_r_interface
# #' @param .data A dataset containing \code{time} and \code{conc}.  Can
# #'   also be a \code{pulse_sim} object.
# #' @param time A string. Name of the time variable in \code{data}.
# #' @param conc A string. Name of the hormone concentration variable in
# #'   \code{data}.
# #' @param spec An object of class \code{spec}, created by
# #'   \code{pulse_spec()}, specifying the priors, starting values, and proposal
# #'   variances to use.
# #' @param iters Number of iterations for MCMC
# #' @param thin Thinning to apply to MCMC chains (i.e. Keep every 'thin'th
# #'   sample).
# #' # param burnin Burn-in to apply to MCMC chains (i.e. remove first 'burnin'
# #' # samples). Applied prior to thinning.
# #' @param use_tibble Return chains as tbl_df class data frames, from the tibble
# #'   package.  Mostly used for the print.tbl_df method, which limits the rows and
# #'   columns printed to those which fit in the console.
# #' @param verbose Prints diagnostics and estimates every 5000th iteration.
# #'   Default is \code{FALSE}.
# #' @import tibble
# #' @keywords pulse fit
# #' @examples
# #' this_pulse <- simulate_pulse()
# #' this_spec  <- pulse_spec()
# #' this_fit   <- fit_pulse(.data = this_pulse, iters = 1000, thin = 10,
# #'                         spec = this_spec)
# #' @export





re_chain <- function(fit) {
  # extract pulse-specific parameter chain
}
fe_chain <- function(fit) {
  # extract common parameters chain
}

