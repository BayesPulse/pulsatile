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
install("../pulsatile")

library(pulsatile)

# Read in dataset
dat <- read_delim("test/pulse_reference_001.dat", delim = " ") %>% tbl_df %>%
  select(-observation)

model_spec <- pulse_spec(.data = dat, iterations = 250000)
model_spec

# Test C functions
#show_args(.data = dat)

#
# ---- Simple test for repeatable results ---- 
#
start_time <- proc.time()
set.seed(999999)
fit_round1 <- fit_pulse(model_spec)
stop_time <- proc.time()
time_round1 <- (stop_time - start_time)/60

start_time <- proc.time()
set.seed(999999)
fit_round2 <- fit_pulse(model_spec)
stop_time <- proc.time()
time_round2 <- (stop_time - start_time)/60

all(fit_round1[[1]] == fit_round2[[1]])
time_round1 == time_round2





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


