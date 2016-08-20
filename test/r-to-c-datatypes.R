#-------------------------------------------------------------------------------
# Testing getting data to the C functions
#
#-------------------------------------------------------------------------------
options(scipen = 99)
setwd("~/Projects/Rpackages/pulsatile/")
library(readr)
library(dplyr)
library(pryr)
library(devtools)
library(roxygen2)

document()
install("../pulsatile")

library(pulsatile)

# Read in dataset
dat <- read_delim("test/pulse_reference_001.dat", delim = " ") %>% tbl_df %>%
  select(-observation)

model_spec <- pulse_spec(.data = dat, iterations = 5000)
model_spec

# Test C functions
#test_inout(x = dat)
#test_inout(x = model_spec)
#show_args(.data = dat)
set.seed(999999)
fit <- fit_pulse(model_spec)

#dat %>% str

pulses <- fit[[2]] %>% do.call(rbind, .) %>% as.data.frame %>% tbl_df  

