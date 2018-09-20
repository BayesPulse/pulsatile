
library(tidyverse)
library(devtools)
document()
# check()
install("./")


source("../libpulsatile/R-package/data-raw/test_data.R")
library(pulsatile)


set.seed(2018-09-18)
spec <- pulsatile::pulse_spec(location_prior_type = "strauss", #order-statistic") #, 
                              prior_location_gamma = 0,
                              prior_location_range = 40)
fit <- pulsatile::fit_pulse(test_data, spec = spec, iters = 250000, thin = 50, burnin = 100000,
                            verbose = TRUE)
pulsatile::bp_location_posterior(fit) + ggtitle("PULSATILE")
pulsatile::bp_trace(fit) + ggtitle("PULSATILE")
pulsatile::pulse_chain(fit) %>%
  group_by(iteration) %>%
  mutate(timediff = location - lag(location)) #%>%
#   filter(pulse_num == 3) 

plot(test_data) +
  geom_vline(data = test_data$parameters, aes(xintercept = location)) +
  ggtitle("PULSATILE") 


bp_trace(fit)

