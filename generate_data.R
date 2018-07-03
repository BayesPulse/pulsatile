#-------------------------------------------------------------------------------
# generate_data.R
#   Generate simulated data for debugging single subject model
#
# Author: M.Mulvahill
# Notes:
#   - generating 4 series, 3 clean and 1 noisy
# 
#   - Other than changing mass_mean to 1.25 and mass_sd to 0.5 in the noisy
#   series, the values for simulation are defaults:
#   baseline              = 2.6
#   mean pulse mass       = 3.5
#   SD/var of pulse mass  = 1.6 (sd)
#   mean pulse width      = 5
#   SD/var of pulse width = 1 (sd)
#   half-life             = 45
#   model error           = 0.005
# 
#-------------------------------------------------------------------------------

setwd("~/Projects/BayesPulse/Software/singlesubject_debug")

# Install/load packages
if (!require(devtools))  install.packages("devtools")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(pulsatile)) devtools::install_github("bayespulse/pulsatile")
library(tidyverse)
library(pulsatile)

# Simulate series
set.seed(2018-05-02)

three_clean_series <- lapply(1:3, function(x) simulate_pulse())
one_noisy_series   <- simulate_pulse(mass_mean = 1.25, mass_sd = 0.5)

plot(one_noisy_series)


# Save datasets and true parameters


lapply(1:3, function(x) {
  write_delim(three_clean_series[[x]]$data, #%>% select(time, concentration), 
              path = paste0("./data/clean_series_", x, ".dat"),
              delim = " ", col_names = FALSE)
})
lapply(1:3, function(x) {
  write_delim(three_clean_series[[x]]$parameters, 
              path = paste0("./data/clean_series_", x, "_trueparameters.dat"),
              delim = " ", col_names = FALSE)
})


lapply(1, function(x) {
  write_delim(one_noisy_series$data, #%>% select(time, concentration), 
              path = paste0("./data/noisy_series_", x, ".dat"),
              delim = " ", col_names = FALSE)
})
lapply(1, function(x) {
  write_delim(one_noisy_series$parameters, 
              path = paste0("./data/noisy_series_", x, "_trueparameters.dat"),
              delim = " ", col_names = FALSE)
})


# Save R versions for easier plotting and viewing of true parameters
clean_series_1 <- three_clean_series[[1]]
clean_series_2 <- three_clean_series[[2]]
clean_series_3 <- three_clean_series[[3]]
rm(three_clean_series)
save.image("./data/sim_series_r_versions.RData")


# TODO:
# - copy problem patient data here (p24.l?) and run with code and write script for
#   viewing results (is pulse locaiton really that bad?)

