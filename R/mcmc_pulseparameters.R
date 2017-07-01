library(pulsatile)
library(purrr)

# 
as_pulse <- function(location = 100, 
                     mass = rnorm(1, 3.5, 1),
                     width = rnorm(1, 35, 1), 
                     mass_kappa = rnorm(1, 1, 0.25),
                     width_kappa = rnorm(1, 1, 0.25),
                     contribution_index = NULL) {

  pulse_parms <- c("location"  = location, 
                  "mass"      = mass,
                  "width"     = width,
                  "mass_kappa"  = mass_kappa,
                  "width_kappa" = width_kappa)

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
  sd_mass    <- 1.5
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





mymcmc <- function(iters = 10000) {

  this_pulse  <- simulate_pulse(num_obs = 144,
                                interval = 10,
                                mass_sd = 1.5)
  test_pulses <- this_pulse$parameters %>% 
    pmap(list) %>% map(~ as_pulse(location    = .x$location,
                                  mass        = .x$mass,
                                  width       = .x$width,
                                  mass_kappa  = .x$mass_kappa,
                                  width_kappa = .x$width_kappa)) %>%
    do.call(as_pulse_list, .)
#   test_pulses <- as_pulse_list(as_pulse(), as_pulse(), as_pulse()) 
  proposal_sd <- 20
  priors <- pulse_spec(prior_max_sd_mass = 5) %>% .$priors
  rtn_sdmass <- NULL
  accepted <- 0

  for (i in 1:iters) {

    rtn_sdmass <- c(rtn_sdmass, test_pulses$parms$sd_mass)
    result     <- cdraw_re_sd(pulse_list = test_pulses,
                              priors = priors,
                              proposal_sd = proposal_sd)

    test_pulses$parms$sd_mass <- result[["value"]]
    accepted    <- (accepted + result[["accepted"]])
    accept_rate <- accepted/i

    
  }

  cat("Final Acceptance Rate:", accept_rate, "\n")
  return(list("df" = data.frame(i = 1:iters, 
                                sdmass = rtn_sdmass, 
                                accept_rate = accept_rate),
              "sim" = this_pulse))

}


################################################################################
# Draw standard deviation of the random effects -- currently for mass only
################################################################################
draw_re_sd <- function(pulse_list,
                       priors,
                       proposal_sd) {

  current_value <- pulse_list$parms$sd_mass
  prior_max <- priors$max_sd$mass

  # draw new value
  new_value <- rnorm(1, current_value, proposal_sd)

  if (new_value < 0 | new_value > prior_max) {
    return(list("value" = current_value, "accepted" = 0))
  }

  parms <- pulse_list$parms
  pulses <- map(pulse_list$pulses, ~ .x$parameters) %>% purrr::transpose(.) %>%
    map(~ do.call(c, .x))

  first_part  <- pulse_list$num_pulses * (log(current_value) - log(new_value))
  second_part <- pulse_list$num_pulses * ((1 / current_value^2) - (1 / new_value^2))
  third_part  <- sum(0.5 *  (pulses$mass - parms$mean_mass) ^ 2 * pulses$mass_kappa)
  old_int     <- sum(pnorm(pulse_list$mean_mass * sqrt(pulses$mass_kappa) / current_value, log.p = TRUE))
  new_int     <- sum(pnorm(pulse_list$mean_mass * sqrt(pulses$mass_kappa) / new_value, log.p = TRUE))

  rho <- first_part + second_part * third_part - new_int + old_int
  rho <- min(0, rho)

  if (log(runif(1, min = 0, max = 1)) < rho) {
    return(list("value" = new_value, "accepted" = 1))
  } else {
    return(list("value" = new_value, "accepted" = 0))
  }

}

library(parallel)
library(compiler)
cdraw_re_sd <- cmpfun(draw_re_sd)

# set.seed(1)
# test1 <- mymcmc(10000)
# test1 %>% filter((i %% 50) == 1) %>% ggplot(aes(x = i, y = sdmass)) + geom_path()
# mean(test1$sdmass)
# sd(test1$sdmass)

set.seed(999)
test2 <- mclapply(1:100, function(.x) { 
                    c("run" = .x, mymcmc(10000))
                       }, mc.cores = detectCores())
test2_dat <- test2 %>% map(~ cbind("run" = .x$run, .x$df)) %>% do.call(rbind, .) %>% as_data_frame
test2_dat <- test2_dat %>% filter((i%%10) == 1)
test2_dat %>% group_by(run) %>% summarise(mean = mean(sdmass))
test2_dat %>% summarise(mean = mean(sdmass))








