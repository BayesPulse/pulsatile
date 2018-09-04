#-------------------------------------------------------------------------------
# Debugging workflow for Nichole
#
#-------------------------------------------------------------------------------


library(pulsatile)
library(purrr)
library(furrr)


#---------------------------------------
# Simulate dataset
#---------------------------------------
set.seed(999)

# IPI of 5
simulate_data <- function(sampling_interval) {

  simulate_pulse(num_obs = 144,
                 interval = sampling_interval,
                 error_var = 0.005,
                 ipi_mean = 12,
                 ipi_var = 40,
                 ipi_min = 4,
                 mass_mean = 3.5,
                 mass_sd = 1.6,
                 width_mean = 5,
                 width_sd = 1,
                 constant_halflife = 45,
                 constant_baseline = 2.6,
                 halflife_mean = NULL,
                 halflife_var = NULL,
                 baseline_mean = NULL,
                 baseline_var = NULL)
}

# ?simulate_pulse # view help file for arguments
sims_ipi10 <- map(1:20, function(x) { simulate_data(sampling_interval = 10) })
sims_ipi5  <- map(1:20, function(x) { simulate_data(sampling_interval = 5) })
sims_ipi1  <- map(1:20, function(x) { simulate_data(sampling_interval = 1) })


#---------------------------------------
# Create model specification used in all fits
#---------------------------------------
myspec <- pulse_spec(location_prior_type = "order-statistic", # or "strauss"
                     prior_mass_mean = 3.5,
                     prior_mass_var = 100,
                     prior_width_mean = 5,
                     prior_width_var = 100,
                     prior_baseline_mean = 2.6,
                     prior_baseline_var = 100,
                     prior_halflife_mean = 45,
                     prior_halflife_var = 100,
                     prior_error_alpha = 1e-04,
                     prior_error_beta = 1e-04,
                     prior_location_gamma = NULL,
                     prior_location_range = NULL,
                     prior_max_sd_mass = 100,
                     prior_max_sd_width = 150,
                     prior_mean_pulse_count = 12,
                     sv_mass_mean = 3.5,
                     sv_width_mean = 5,
                     sv_baseline_mean = 2.6,
                     sv_halflife_mean = 45,
                     sv_error_var = 0.005,
                     sv_mass_sd = 1.6,
                     sv_width_sd = 1)


#---------------------------------------
# Fit models
#---------------------------------------
plan(multiprocess)

list_of_sims <- c(sims_ipi10, sims_ipi5, sims_ipi1)
list_of_results <-
  future_map(list_of_sims,
             ~ fit_pulse(.x, time = "time", conc = "concentration",
                         spec = myspec, iters = 250000, thin = 50, burnin = 100000,
                         verbose = TRUE))

# SAVED SIM DATA AND CHAINS IN ONE OF OUR SHARED DROPBOX FOLDERS
save(list_of_sims, list_of_results, 
     file = "~/Dropbox/Work/HormoneProjects/Hormone_Code_Versionn/simulated_data_for_nichole_withfix.RData")




