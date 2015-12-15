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

test_inout(as.matrix(dat))
test_inout(dat)

show_args(dat)

dat %>% str


