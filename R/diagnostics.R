#-------------------------------------------------------------------------------
# diagnostics.R
#     - Diagnostic plots and summary stats
#-------------------------------------------------------------------------------

#--------------------------------------------
# STAN examples
#   Some examples http://mc-stan.org/bayesplot/
#--------------------------------------------
# library(bayesplot)
# library(rstanarm)
# library(ggplot2)
# fit        <- stan_glm(mpg ~ ., data = mtcars)
# posterior  <- as.matrix(fit)
# 
# plot_title <- ggtitle("Posterior distributions with medians and 80% intervals")
# mcmc_areas(posterior, pars = c("cyl", "drat", "am", "wt"), prob = 0.8) +
#   plot_title
# ppc_intervals(y = mtcars$mpg, yrep = posterior_predict(fit), x = mtcars$wt, prob = 0.5) +
#   labs(x = "Weight (1000 lbs)",
#        y = "MPG",
#        title = "50% posterior predictive intervals \nvs observed miles per gallon",
#        subtitle = "by vehicle weight") +
#   panel_bg(fill = "gray95", color = NA) +
#   grid_lines(color = "white")

#--------------------------------------------
# STAN and other Bayesian R package functions to implement
#
# Posterior predicted values/plot
#   - rstanarm::posterior_predict()
#   - rstanarm::ppc_dens_overlay()
#   - rstanarm::ppc_intervals()
# Posterior densities
#   - rstanarm::mcmc_areas()
#

#--------------------------------------------
# Other functions w/ no corollary from major package
#
#   Extract heirarchical chains
#     - hormone1_chain() ??? driver_chain()
#     - hormone2_chain() ??? response_chain()
#     - pulse_chain()
#     - subject_chain()
#     - pop_chain()
#     - chains() --> list/str chains available and functions to extract them
#
#   Summary -- follow rstanarm::summary.stanreg's example


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


# mcmc_trace <- function() {}
# mcmc_posteriors <- function() {}
# mcmc_locations <- function() {}
# 
# 
# 
#   temp <- 
#     common %>%
#       gather(key = key, value = value, num.pulses:sd.widths) %>%
#       filter(dataset %in% dataset.nums) %>%
#       group_by(dataset) %>%
#       do( 
#         # Trace plots
#         trace.figs = 
#         {
#           trace.fig <- 
#             ggplot(., aes(x = iteration, y = value)) +
#               geom_path(size = 0.10) +
#               facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free") +
#               ggtitle(paste("Dataset", ifelse(is.null(.$orig.dataset), 
#                                               unique(.$dataset),
#                                               unique(.$orig.dataset))))
#           #print(trace.fig)
#         },
#         # Posterior densities
#         post.figs = 
#         {
#           post.fig <- 
#             ggplot(., aes(x = value)) +
#               geom_histogram(aes(y = ..density..), size = 0.15) +
#               #geom_density(alpha=.2, fill="#FF6666") +
#               facet_wrap( ~ key, ncol = 2, nrow = 4, scales = "free") +
#               ggtitle(paste("Dataset", ifelse(is.null(.$orig.dataset), 
#                                               unique(.$dataset),
#                                               unique(.$orig.dataset))))
#           #suppressMessages(print(post.fig))
#         }
#       )
# 
#   location.figs.lst <- 
#     pulse %>%
#       filter(dataset %in% dataset.nums) %>%
#       group_by(dataset) %>%
#       do(
#         # Location histograms
#         location.figs = 
#         {
#           location.fig <- 
#             ggplot(., aes(x = location)) +
#               geom_histogram(binwidth = 5) +
#               theme(panel.grid.minor = element_line(colour="lightgrey", size=0.5)) + 
#               theme(panel.grid.major = element_line(colour="lightgrey", size=0.5)) + 
#               scale_x_continuous(breaks = seq(-50, max.time+50, 50),
#                                  minor_breaks = seq(-50, max.time+50, 10),
#                                  limits = c(-50, max.time+50)) 
#         }
#       )
# 
#   sim.figs.lst <-
#     sim %>%
#     filter(dataset %in% dataset.nums) %>%
#     group_by(dataset) %>%
#     do(
#        sim.figs = 
#        {
#          sim.fig <-
#            ggplot(., aes(x = time, y = concentration)) +
#              geom_path() +
#              geom_point() + 
#              theme(panel.grid.minor = element_line(colour="lightgrey", size=0.5)) + 
#              theme(panel.grid.major = element_line(colour="lightgrey", size=0.5)) + 
#              scale_x_continuous(breaks = seq(-50, max.time+50, 50),
#                                 minor_breaks = seq(-50, max.time+50, 10),
#                                 limits = c(-50, max.time+50)) 
#        }
#     )
