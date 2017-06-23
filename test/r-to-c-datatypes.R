#-------------------------------------------------------------------------------
# Testing getting data to the C functions
#
#-------------------------------------------------------------------------------
#options(scipen = 99)
setwd("~/Projects/BayesPulse/pulsatile/")

# library(dplyr)
# library(tidyr)
#library(pryr)
library(devtools)
# library(roxygen2)
library(ggplot2)
library(magrittr)
library(ggthemes)
theme_set(theme_tufte())

devtools::document()
devtools::check()
devtools::install("../pulsatile", build_vignettes = TRUE)

library(pulsatile)

# NOTE: eta is not working right -- often -nan or huuuuuge
set.seed(9999)
this_pulse <- simulate_pulse()
model_spec <- pulse_spec(location_prior_type = "order-statistic")
fit_test   <- fit_pulse(.data = this_pulse, iters = 25000, thin = 50,
                        spec = model_spec, verbose = TRUE)
str(fit_test)

plot(this_pulse)

?chains
chains(fit_test)
pulse_chain(fit_test)
common_chain(fit_test)

fit_test$common_chain %>% 
  ggplot(aes(x = iteration, y = mean_pulse_width)) + geom_path()
# will add summary and print s3 methods
#summary(this_pulse)
this_pulse

traceplots <- 
  fit_test$common_chain %>% 
  gather(key = parameter, value = value, -iteration) %>%
  ggplot(aes(x = iteration, y = value)) +
    geom_path() +
    facet_wrap(~ parameter, ncol = 3, scales = "free")

posterior_dens <- 
  fit_test$common_chain %>% 
  gather(key = parameter, value = value, -iteration) %>%
  ggplot(aes(x = value)) +
    geom_histogram() +
    facet_wrap(~ parameter, ncol = 3, scales = "free")


##############################
# Various prior parameter sets to fit model with
##############################
model_spec <- pulse_spec(location_prior_type = "order-statistic")
model_spec_strauss <- pulse_spec(location_prior_type  = "strauss",
                                 prior_location_range = 40,
                                 prior_location_gamma = 0)


#
# ---- Simple test for repeatable results ---- 
#
n_iters <- 100000
n_thin  <- 50

start_time <- proc.time()
set.seed(999999)
fit_round1 <- fit_pulse(.data = this_pulse,
                        iters = n_iters,
                        thin  = n_thin,
                        spec  = model_spec,
                        verbose = TRUE)
stop_time <- proc.time()
time_round1 <- (stop_time - start_time)/60

start_time <- proc.time()
set.seed(999999)
fit_round2 <- fit_pulse(.data = this_pulse,
                        iters = n_iters,
                        thin  = n_thin,
                        spec  = model_spec,
                        verbose = TRUE)
stop_time <- proc.time()
time_round2 <- (stop_time - start_time)/60

start_time <- proc.time()
set.seed(999999)
fit_strauss <- fit_pulse(.data = this_pulse,
                         iters = n_iters,
                         thin  = n_thin,
                         spec  = model_spec_strauss,
                        verbose = TRUE)
stop_time <- proc.time()
time_round3 <- (stop_time - start_time)/60



########################################
# Compare results
########################################
all(fit_round1[["common_chain"]] == fit_round2[["common_chain"]])
all(fit_round1[["pulse_chain"]] == fit_round2[["pulse_chain"]])
time_round1 == time_round2


########################################
# Sanity check on pulse count
########################################
hist(fit_round1$common_chain$num_pulses) 
hist(fit_round2$common_chain$num_pulses) 
hist(fit_strauss$common_chain$num_pulses) 

dev.off()
dev.off()
dev.off()

########################################
# Check Strauss 
########################################
library(dplyr)
library(tidyr)
library(ggplot2)

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
    geom_path(data = fit_round1$data$data, aes(x = time, y = concentration)) 

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





# check new vs old simulation function
# source("R/simulate.R")
# source("../simulate.R")
# set.seed(999)
# after <- new_simulate_pulse()
# set.seed(999)
# before <- simulate_pulse()
# 
# identical(new_simulate_pulse, simulate_pulse)
# identical(after$pulse_data, before$pulse_data)
# identical(after$pulse_parms, before$pulse_parms)
# dplyr::full_join(after$pulse_data, before$pulse_data, by = c("observation", "time"))
# cbind(after$pulse_parms, before$pulse_parms)
