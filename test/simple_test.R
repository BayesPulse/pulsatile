setwd("~/Projects/Rpackages/pulsatile")
library(readr)
library(dplyr)
library(pulsatile)

# Read in dataset
dat <- read_delim("test/pulse_reference_001.dat", delim = " ") %>% tbl_df %>%
  select(-observation)
save(dat, file = "test/data.Rdata")

model_spec <- pulse_spec(.data = dat)
model_spec

# Test C functions
fit_pulse(model_spec)



