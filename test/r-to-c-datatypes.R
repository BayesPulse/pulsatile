#-------------------------------------------------------------------------------
# Testing getting data to the C functions
#
#-------------------------------------------------------------------------------

library(pulsatile)

library(readr)
library(dplyr)

# Read in dataset
dat <- read_delim("test/pulse_reference_001.dat", delim = " ") %>% tbl_df


