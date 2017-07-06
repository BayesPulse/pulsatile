# 
# #
# #
# #
# #
# 
# 
as_pulse <- function(location  = 100, 
                     mass      = rnorm(1, 3.5, 1),
                     width     = rnorm(1, 35, 10), 
                     mass_eta  = 1, #rnorm(1, 0.5, 0.25),
                     width_eta = 1, #rnorm(1, 0.5, 0.25),
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
  sd_width   <- 10
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



mymcmc <- function() {


  test_pulses <- as_pulse_list(as_pulse(), as_pulse(), as_pulse(), as_pulse(), 
                               as_pulse(), as_pulse(), as_pulse(), as_pulse(),
                               as_pulse(), as_pulse(), as_pulse(), as_pulse(),
                               as_pulse(), as_pulse(), as_pulse(), as_pulse(),
                               as_pulse(), as_pulse(), as_pulse(), as_pulse())
  proposal_sd <- 5
  priors <- pulse_spec(prior_max_sd_mass = 10) %>% .$priors
  rtn_sdmass <- 0

  for (i in 1:25000) {

    rtn_sdmass <- c(rtn_sdmass, test_pulses$parms$sd_mass)
    result <- draw_sd_between_pulses(pulse_list = test_pulses,
                                     priors = priors,
                                     proposal_sd = proposal_sd)

    test_pulses$parms$sd_mass <- result$value
  }

  return(rtn_sdmass)


}

# prior max = 10
# proposal sd = 15
# sd mass = 1
#set.seed(1)
test <- mymcmc()
test <- data_frame("i" = 1:length(test), "sd" = test) 
ggplot(test, aes(x = i, y = sd)) + geom_path()
# mean(sqrt((test$sd^2)/4))
mean(test$sd)
median(test$sd)



