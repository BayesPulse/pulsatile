#-------------------------------------------------------------------------------
# Testing getting data to the C functions
#
#-------------------------------------------------------------------------------
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

model_spec <- pulse_spec()
model_spec

test_inout(x = as.matrix(dat)) #, specification = model_spec)
test_inout(x = dat)

test_inout(x = model_spec)


show_args(.data = dat)

dat %>% str


