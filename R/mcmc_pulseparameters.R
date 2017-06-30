# 
# #
# #
# #
# #
# 
# 
as_pulse <- function(location = 100, 
                     mass = rnorm(1, 3.5, 1),
                     width = rnorm(1, 35, 1), 
                     mass_eta = rnorm(1, 1, 0.25),
                     width_eta = rnorm(1, 1, 0.25),
                     contribution_index = NULL) {

  pulse_parms <- c("location"  = location, 
                  "mass"      = mass,
                  "width"     = width,
                  "mass_eta"  = mass_eta,
                  "width_eta" = width_eta)

  structure(list("parameters" = pulse_parms,
                 "contribution_index" = contribution_index),
            class = c("pulse"))

}

#' as_pulse_list
#' 
#' @param ... bunch o' pulses
as_pulse_list <- function(..., .data = NULL, num_obs = 144, sampling_interval = 10) {

  pulses <- list(...) # check all of pulse class
  num_pulses <- length(pulses)
  duration <- num_obs * sampling_interval # total minutes

  baseline_concentration <- 2.6
  half_life  <- 45
  error_sd   <- 0.25
  mean_mass  <- 3.5
  mean_width <- 35
  sd_mass    <- 1
  sd_width   <- 1
  decay      <- log(2) / half_life
  parms <- list("baseline_concentration" = baseline_concentration, 
                "half_life"  = half_life,
                "mean_mass"  = mean_mass,
                "mean_width" = mean_width,
                "sd_mass"    = sd_mass,
                "sd_width"   = sd_width,
                "decay"      = decay)

#   mean_concentration
  structure(list("parms" = parms, "num_pulses" = num_pulses, 
                 "pulses" = pulses),
            class = c("pulse_list"))
}




draw_re_sd <- function(pulse_list,
                       priors,
                       proposal_sd) {

  current_value <- pulse_list$parms$sd_mass
  prior_max <- priors$max_sd$mass

  # draw new value
  new_value <- rnorm(1, current_value, proposal_sd)

  if (new_value < 0 | new_value > prior_max) {
    return(current_value)
  }

  pulses <- map(pulse_list$pulses, ~ .x$parameters) %>% purrr::transpose(.) %>%
    map(~ do.call(c, .x))

  first_part <- 
    pulse_list$num_pulses * (log(current_value) - log(new_value))
  second_part <- 
    pulse_list$num_pulses * ((1 / current_value^2) - (1 / new_value^2))
  third_part <-
    sum(0.5 *  (pulses$mass - parms$mean_mass) ^ 2 * pulses$mass_eta)
  old_int <- 
    sum(pnorm(pulse_list$mean_mass * sqrt(pulses$mass_eta) / current_value, log.p = TRUE))
  new_int <- 
    sum(pnorm(pulse_list$mean_mass * sqrt(pulses$mass_eta) / new_value, log.p = TRUE))

  rho <- first_part + second_part * third_part - new_int + old_int
  rho <- min(0, rho)

  if (log(runif(1, min = 0, max = 1)) < rho) {
    return(new_value)
  } else {
    return(current_value)
  }

}

mymcmc <- function() {


  test_pulses <- as_pulse_list(as_pulse(), as_pulse(), as_pulse()) 
  proposal_sd <- 15
  priors <- pulse_spec(prior_max_sd_mass = 10) %>% .$priors
  rtn_sdmass <- 0

  for (i in 1:100000) {

    rtn_sdmass <- c(rtn_sdmass, test_pulses$parms$sd_mass)
    test_pulses$parms$sd_mass <- draw_re_sd(pulse_list = test_pulses,
                                            priors = priors,
                                            proposal_sd = proposal_sd)
    
  }

  return(rtn_sdmass)


}




