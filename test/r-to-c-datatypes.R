#-------------------------------------------------------------------------------
# Testing getting data to the C functions
#
#-------------------------------------------------------------------------------
options(scipen = 99)
setwd("~/Projects/Rpackages/pulsatile/")
library(readr)
library(dplyr)
library(tidyr)
library(pryr)
library(devtools)
library(roxygen2)
library(ggplot2)
library(ggthemes)
theme_set(theme_tufte())

document()
install("../pulsatile", build_vignettes = TRUE)

library(pulsatile)
#data(simpulse_reference001)

sim <- simulate_pulsets()
dat <- sim$pulse_data

# Replace by above^
# Read in dataset
#dat <- read_delim("test/pulse_reference_001.dat", delim = " ") %>% 
#  tbl_df %>%
#  select(-observation)

model_spec <- pulse_spec(location_prior_type = "order-statistic")
model_spec
model_spec_strauss <- pulse_spec(location_prior_type = "strauss",
                                 prior_location_range = 40,
                                 prior_location_gamma = 0)
model_spec_strauss


#
# ---- Simple test for repeatable results ---- 
#
n_iters <- 5000
n_thin  <- 1

start_time <- proc.time()
set.seed(999999)
fit_round1 <- fit_pulse(data           = dat,
                        iterations     = n_iters,
                        thin           = n_thin,
                        pulse_spec_obj = model_spec)
stop_time <- proc.time()
time_round1 <- (stop_time - start_time)/60

start_time <- proc.time()
set.seed(999999)
fit_round2 <- fit_pulse(data           = dat,
                        iterations     = n_iters,
                        thin           = n_thin,
                        pulse_spec_obj = model_spec)
stop_time <- proc.time()
time_round2 <- (stop_time - start_time)/60

start_time <- proc.time()
set.seed(999999)
fit_strauss <- fit_pulse(data           = dat,
                         iterations     = n_iters,
                         thin           = n_thin,
                         pulse_spec_obj = model_spec_strauss)
stop_time <- proc.time()
time_round3 <- (stop_time - start_time)/60



########################################
# Compare results
########################################
all(fit_round1[[3]] == fit_round2[[3]])
all(fit_round1[[4]] == fit_round2[[4]])
time_round1 == time_round2


########################################
# Sanity check on pulse count
########################################
hist(fit_round1$common_chain$num_pulses) 
# NOTE: Somethings up with the birth-death process.  Think it started with the
# addition of the location_prior_type option in pulse_spec() and in C code.

########################################
# Check Strauss 
########################################
all((fit_round1[[3]] == fit_strauss[[3]]))
# min. distance between locations should be approx, but > 40
fit_strauss$pulse_chain %>% 
  group_by(iteration) %>%
  mutate(distance = location - lag(location, k = 1)) %>%
  ungroup %>% with(., min(distance, na.rm=T))



########################################
# Visualize results
########################################

pulses <- fit_round1$pulse_chain # %>% do.call(rbind, .) %>% as.data.frame %>% tbl_df  
common <- fit_round1$common_chain # %>% as.data.frame %>% tbl_df %>% mutate(iteration = 1:n()) %>% select(iteration, everything())

timeseries <- 
  ggplot() +
    geom_path(data = fit_round1$data, aes(x = time, y = concentration)) 

location_hist <- 
  ggplot() + geom_histogram(data = as.data.frame(pulses), aes(x = location, y = ..density..))

traceplots <- 
  common %>% gather(key = parameter, value = value, -iteration) %>%
  ggplot(aes(x = iteration, y = value)) +
    geom_path() +
    facet_wrap(~ parameter, ncol = 3, scales = "free")


########################################
# inspect object size 
########################################
# comparing matrix to data frame to tibble
object.size(pulses)
test1 <- as.data.frame(pulses)
test2 <- as_data_frame(test1)
object.size(test1)
object.size(test2)


