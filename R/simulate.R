
#' Simulate pulsatile hormone data
#' 
#' @description \code{\link{simulate_pulse}} simulates a time series dataset
#'   representing blood concentration measurements of a pulsatile hormone. Both
#'   the time series and a dataset of the individual pulse characteristics are
#'   returned. 
#'
#' @param num_of_observations Number of observations to simulate.  Duration of
#'   observation window equals \code{num_of_observations} times
#'   \code{sampling_interval}.
#' @param sampling_interval Time between observations, typically 6-10 minutes.
#' @param error_var Variance of the error added at each observation, ~ N(0, sqrt(error_var)).
#' @param ipi_mean Mean number of sampling units between pulses (mean inter-pulse interval).
#' @param ipi_var Variance of gamma for drawing interpulse interval
#' @param ipi_min Minimum number of units between pulses
#' @param mass_mean Mean pulse mass
#' @param mass_sd Standard deviation of pulse mass
#' @param width_mean Mean pulse width
#' @param width_sd Standard deviation of pulse width
#' @param halflife_mean Mean of half-life
#' @param halflife_var Variance of half-life
#' @param baseline_mean Mean of baseline
#' @param baseline_var Variance of baseline
#' @param constant_halflife To use a constant (specified) half-life, set this
#'   to a constant [0,inf). Mean and variance of half-life are not used if this
#'   is non-null.
#' @param constant_baseline To use a constant (specified) baseline, set this to
#'   a constant [0,inf). Mean and variance of baseline are not used if this is
#'   non-null.
#' @return A object of class \code{pulse_sim} containing time-series dataset
#'   and dataset of characteristics of each pulse
#' @keywords pulse simulation
#' @examples
#' this_pulse <- simulate_pulse()
#' str(this_pulse)
#' plot(x = this_pulse$pulse_data$time, 
#'      y = this_pulse$pulse_data$concentration,
#'      type = 'l')
#' # Add this -> summary(this_pulse)
#' # Add this -> plot()
#' # Add this -> print()
#' @export
simulate_pulse <- function(num_of_observations = 144,
                           sampling_interval   = 10,
                           error_var           = 0.005,
                           ipi_mean            = 12,
                           ipi_var             = 40,
                           ipi_min             = 4,
                           mass_mean           = 3.5,
                           mass_sd             = 1.6,
                           width_mean          = 35,
                           width_sd            = 5,
                           halflife_mean       = NULL,
                           halflife_var        = NULL,
                           baseline_mean       = NULL,
                           baseline_var        = NULL,
                           constant_halflife   = 45,
                           constant_baseline   = 2.6) {

  #---------------------------------------
  # Helper functions for drawing hormone concentration (w/o error)
  # Normal CDF 
  erfFn <- function(x) 2 * stats::pnorm(x * sqrt(2), 0, 1) - 1

  # Model concentration over time given pulse parameters  
  meanI <- function(sampling_interval, b, a, tau1, lam, s2p){
    b + (a / 2) * 
      exp((tau1 - sampling_interval) * lam + 0.5 * lam^2 * s2p) * 
      (1 + erfFn((sampling_interval - (tau1 + lam * s2p)) / sqrt(2 * s2p)))
  }

  #---------------------------------------
  # Get baseline concentration
  if (is.null(constant_baseline)) {
    while(B <= 0) B <- stats::rnorm(1, baseline_mean, sqrt(baseline_var))
  } else {
    B <- constant_baseline
  }

  #---------------------------------------
  # Get half-life of hormone
  #   H is half-life, H=ln(2)/lambda_x, where lambda_x is the decay constant
  if (is.null(constant_halflife)) {
    while(H <= 8) H <- stats::rnorm(1, halflife_mean, sqrt(halflife_var))
  } else {
    H <- constant_halflife
  }

  #---------------------------------------
  # Generate pulse locations
  #   Using a renewal process, define by interpulsatile interval and variance
  #   then convert gamma parameters
  # - mean = alpha / beta
  # - var = alpha / beta^2
  gammamean <- ipi_mean - ipi_min
  alpha     <- gammamean * gammamean / ipi_var
  beta      <- gammamean / ipi_var

  tau <- rep(0, 25)
  tau[1] <- sampling_interval * (stats::rgamma(1, alpha, beta) - 10)

  i <- 1
  while (tau[i] < (num_of_observations * sampling_interval)){
    i      <- i + 1
    tmp    <- ipi_min + stats::rgamma(1, alpha, beta)
    tau[i] <- tau[i-1] + (sampling_interval * tmp)
  }

  # Reduce pulse location vector to values within time range (<0, 1440)
  tau <- subset(tau, tau < (num_of_observations * sampling_interval))
  tau <- subset(tau, tau != 0)
  np  <- length(tau) # No. of pulses 

  #---------------------------------------
  # Generate pulse-specific parameters 
  A     <- rep(0, np)             # pulse mass
  s2p   <- rep(0, np)             # pulse width
  taxis <- seq(10, (num_of_observations * sampling_interval), sampling_interval)
  ytmp  <- rep(0, length(taxis))  # hormone concentration

  for (i in 1:np) {
    # Log-normal 
    #A[i]   <- exp(stats::rnorm(1, mass_mean, mass_sd))
    #s2p[i] <- exp(stats::rnorm(1, width_mean, width_sd))

    # Truncated T (via gamma normal mixture)
    tvar  <- mass_sd^2  / stats::rgamma(1, 2, 2)        
    t2var <- width_sd^2 / stats::rgamma(1, 2, 2)
    while (A[i] < 0.25)   A[i]   <- stats::rnorm(1, mass_mean, sqrt(tvar))
    while (s2p[i] < 0.5)  s2p[i] <- stats::rnorm(1, width_mean, sqrt(t2var))
  }

  #---------------------------------------
  # Draw mean concentration
  ytmp <- 0
  for (i in 1:np) {
    ytmp   <- ytmp + meanI(taxis, 0, A[i], tau[i], log(2) / H, s2p[i])
  }
  # Add baseline and error
  ysim <- ytmp + B
  errors    <- stats::rnorm((length(taxis)), 0, sqrt(error_var))
  ysimerror <- ysim * exp(errors)

  #---------------------------------------
  # Combine into final simulated datasets
  allpulseparms <- cbind("pulse_no" = seq(1, np),
                         "mass"     = A,
                         "width"    = s2p,
                         "location" = tau)
  allpulseparms <- tibble::as_data_frame(allpulseparms)

  ysim_df   <- cbind("observation"   = 1:length(taxis),
                     "time"          = taxis,
                     "concentration" = ysimerror)
  ysim_df   <- tibble::as_data_frame(ysim_df)

  #---------------------------------------
  # Create return object
  rtn <- structure(list("pulse_data"  = ysim_df, 
                        "pulse_parms" = allpulseparms),
                   class = "pulse_sim")


  return(rtn)

}




#-------------------------------------------------------------------------------
# End of file  
#-------------------------------------------------------------------------------
