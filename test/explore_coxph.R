#-------------------------------------------------------------------------------
# Exploring survival::coxph for handling of call, envir, and formula
#-------------------------------------------------------------------------------

library(survival)
library(dplyr)
library(tibble)

aml <- aml %>% tbl_df

debug(coxph)
fit <- coxph(Surv(time, status) ~ x, data = aml)
undebug(coxph)

