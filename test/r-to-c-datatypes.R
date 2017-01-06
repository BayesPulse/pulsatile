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
data(simpulse_reference001)

# Replace by above^
# Read in dataset
#dat <- read_delim("test/pulse_reference_001.dat", delim = " ") %>% 
#  tbl_df %>%
#  select(-observation)

model_spec <- pulse_spec()
model_spec

# Test C functions
#show_args(.data = dat)

#
# ---- Simple test for repeatable results ---- 
#
start_time <- proc.time()
set.seed(999999)
fit_round1 <- fit_pulse(data           = simpulse_reference001,
                        iterations     = 10000,
                        thin           = 1,
                        pulse_spec_obj = model_spec)
stop_time <- proc.time()
time_round1 <- (stop_time - start_time)/60

start_time <- proc.time()
set.seed(999999)
fit_round2 <- fit_pulse(data = dat, iterations = 10000, thin = 1,
                        pulse_spec_obj = model_spec)
stop_time <- proc.time()
time_round2 <- (stop_time - start_time)/60

all(fit_round1[[3]] == fit_round2[[3]])#[1:1000, ]
all(fit_round1[[4]] == fit_round2[[4]])#[1:1000, ]
time_round1 == time_round2

all((fit_round1[[3]] == fit_round2[[3]])[491:512, ])
#(fit_round1[[4]][[502]] == fit_round2[[4]][[502]])






pulses <- fit[[2]] %>% do.call(rbind, .) %>% as.data.frame %>% tbl_df  
common <- fit[[1]] %>% as.data.frame %>% tbl_df %>% mutate(iteration = 1:n()) %>% select(iteration, everything())

timeseries <- 
  ggplot() +
    geom_path(data = dat, aes(x = time, y = concentration)) +

location_hist <- 
  ggplot() +
  geom_histogram(data = pulses, aes(x = location, y = ..density..))

traceplots <- 
  common %>% gather(key = parameter, value = value, -iteration) %>%
  ggplot(aes(x = iteration, y = value)) +
    geom_path() +
    facet_wrap(~ parameter, ncol = 3, scales = "free")


