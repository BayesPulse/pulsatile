% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit_pulse}
\alias{fit_pulse}
\title{fit_pulse}
\usage{
fit_pulse(.data, time = "time", conc = "concentration", spec,
  iters = 250000, thin = 50, burnin = as.integer(0.1 * iters),
  use_tibble = TRUE, verbose = FALSE)
}
\arguments{
\item{.data}{A dataset containing \code{time} and \code{conc}.  Can
also be a \code{pulse_sim} object.}

\item{time}{A string. Name of the time variable in \code{data}.}

\item{conc}{A string. Name of the hormone concentration variable in
\code{data}.}

\item{spec}{An object of class \code{spec}, created by
\code{pulse_spec()}, specifying the priors, starting values, and proposal
variances to use.}

\item{iters}{Number of iterations for MCMC}

\item{thin}{Thinning to apply to MCMC chains (i.e. Keep every 'thin'th
sample).}

\item{burnin}{Burn-in to apply to MCMC chains (i.e. remove first 'burnin'
samples). Applied prior to thinning.}

\item{use_tibble}{Return chains as tbl_df class data frames, from the tibble
package.  Mostly used for the print.tbl_df method, which limits the rows and
columns printed to those which fit in the console.}

\item{verbose}{Prints diagnostics and estimates every 5000th iteration.
Default is \code{FALSE}.}
}
\description{
Primary function for fitting deconvolution model for time series of pulsatile
hormone data
}
\examples{
this_pulse <- simulate_pulse()
this_spec  <- pulse_spec()
this_fit   <- fit_pulse(.data = this_pulse, iters = 1000, thin = 10, 
                        burnin = 100, spec = this_spec)
}
\keyword{fit}
\keyword{pulse}
