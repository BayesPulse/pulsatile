#-------------------------------------------------------------------------------
# Debugging workflow for Nichole
#
#-------------------------------------------------------------------------------

# Update w/ your path
setwd("~/Repositories/pulsatile/")

# Library for building/compiling/testing packages
library(devtools)

# Other key devtools functions
# devtools::document() # rebuild documentation (and compile source code)
# devtools::check()    # check that package meets CRAN specifications


#---------------------------------------
# Re-compile and install the package
#---------------------------------------
devtools::install("../pulsatile", build_vignettes = FALSE)

#---------------------------------------
# Simulate dataset
#---------------------------------------
set.seed(999)

# ?simulate_pulse # view help file for arguments
sim <- simulate_pulse(num_obs = 144,
                      interval = 10,
                      error_var = 0.005,
                      ipi_mean = 12,
                      ipi_var = 40,
                      ipi_min = 4,
                      mass_mean = 3.5,
                      mass_sd = 1.6,
                      width_mean = 5,
                      width_sd = 1,
                      constant_halflife = 45,
                      constant_baseline = 2.6,
                      halflife_mean = NULL,
                      halflife_var = NULL,
                      baseline_mean = NULL,
                      baseline_var = NULL)

#---------------------------------------
# Create model specification
#---------------------------------------
myspec <- pulse_spec(location_prior_type = "order-statistic", # or "strauss"
                     prior_mass_mean = 3.5,
                     prior_mass_var = 100,
                     prior_width_mean = 5,
                     prior_width_var = 100,
                     prior_baseline_mean = 2.6,
                     prior_baseline_var = 100,
                     prior_halflife_mean = 45,
                     prior_halflife_var = 100,
                     prior_error_alpha = 1e-04,
                     prior_error_beta = 1e-04,
                     prior_location_gamma = NULL,
                     prior_location_range = NULL,
                     prior_max_sd_mass = 100,
                     prior_max_sd_width = 150,
                     prior_mean_pulse_count = 12,
                     sv_mass_mean = 3.5,
                     sv_width_mean = 5,
                     sv_baseline_mean = 2.6,
                     sv_halflife_mean = 45,
                     sv_error_var = 0.005,
                     sv_mass_sd = 1.6,
                     sv_width_sd = 1)


#---------------------------------------
# Fit model
#---------------------------------------
clean1<-read.delim("/Users/nichole/Repositories/singlesubject_debug/data/clean_series_1.dat",sep=" ",header=F,col.names = c("Obs","time","concentration"))
fit <- fit_pulse(clean1, time = "time", conc = "concentration",
                 spec = myspec, iters = 50000, thin = 1, burnin = 0,
                 verbose = TRUE)


#---------------------------------------
# Check chains and posteriors
#---------------------------------------

# traceplots
bp_trace(fit)

# common parms posteriors
bp_posteriors(fit)

# pulse location posterior
bp_location_posterior(fit)



#---------------------------------------
# Fit model
#---------------------------------------
clean2<-read.delim("/Users/nichole/Repositories/singlesubject_debug/data/clean_series_2.dat",sep=" ",header=F,col.names = c("Obs","time","concentration"))
fit2 <- fit_pulse(clean2, time = "time", conc = "concentration",
                 spec = myspec, iters = 50000, thin = 1, burnin = 0,
                 verbose = TRUE)


#---------------------------------------
# Check chains and posteriors
#---------------------------------------

# traceplots
bp_trace(fit2)

# common parms posteriors
bp_posteriors(fit2)

# pulse location posterior
bp_location_posterior(fit2)

mean.mass <-tapply(fit2$pulse_chain$mass,fit2$pulse_chain$iteration,mean)
hist(mean.mass)
mean(mean.mass)
var.mass <-tapply(fit2$pulse_chain$mass,fit2$pulse_chain$iteration,var)
sd.mass<-sqrt(var.mass)
mean(sd.mass)
mean.width <-tapply(fit2$pulse_chain$width,fit2$pulse_chain$iteration,mean)
var.width <-tapply(fit2$pulse_chain$width,fit2$pulse_chain$iteration,var)
mean(mean.width)
sd.width<-sqrt(var.width)
mean(sd.width)


#---------------------------------------
# Create model specification for noisy series
#---------------------------------------
myspec2 <- pulse_spec(location_prior_type = "order-statistic", # or "strauss"
                     prior_mass_mean = 1.25,
                     prior_mass_var = 100,
                     prior_width_mean = 5,
                     prior_width_var = 100,
                     prior_baseline_mean = 2.6,
                     prior_baseline_var = 100,
                     prior_halflife_mean = 45,
                     prior_halflife_var = 100,
                     prior_error_alpha = 1e-04,
                     prior_error_beta = 1e-04,
                     prior_location_gamma = NULL,
                     prior_location_range = NULL,
                     prior_max_sd_mass = 100,
                     prior_max_sd_width = 150,
                     prior_mean_pulse_count = 12,
                     sv_mass_mean = 1.25,
                     sv_width_mean = 5,
                     sv_baseline_mean = 2.6,
                     sv_halflife_mean = 45,
                     sv_error_var = 0.005,
                     sv_mass_sd = 0.5,
                     sv_width_sd = 1)


#---------------------------------------
# Fit model
#---------------------------------------
noisy1<-read.delim("/Users/nichole/Repositories/singlesubject_debug/data/noisy_series_1.dat",sep=" ",header=F,col.names = c("Obs","time","concentration"))
fitn <- fit_pulse(noisy1, time = "time", conc = "concentration",
                 spec = myspec2, iters = 50000, thin = 1, burnin = 0,
                 verbose = TRUE)


#---------------------------------------
# Check chains and posteriors
#---------------------------------------

# traceplots
bp_trace(fitn)

# common parms posteriors
bp_posteriors(fitn)

# pulse location posterior
bp_location_posterior(fitn)

mean.mass <-tapply(fitn$pulse_chain$mass,fitn$pulse_chain$iteration,mean)
hist(mean.mass)
mean(mean.mass)
var.mass <-tapply(fitn$pulse_chain$mass,fitn$pulse_chain$iteration,var)
sd.mass<-sqrt(var.mass)
mean(sd.mass)
mean.width <-tapply(fitn$pulse_chain$width,fitn$pulse_chain$iteration,mean)
var.width <-tapply(fitn$pulse_chain$width,fitn$pulse_chain$iteration,var)
mean(mean.width)
sd.width<-sqrt(var.width)
mean(sd.width)

#---------------------------------------
# Create model specification for P24 real data
#---------------------------------------
myspec24 <- pulse_spec(location_prior_type = "order-statistic", # or "strauss"
                      prior_mass_mean = 1.5,
                      prior_mass_var = 100,
                      prior_width_mean = 5,
                      prior_width_var = 100,
                      prior_baseline_mean = 2.6,
                      prior_baseline_var = 100,
                      prior_halflife_mean = 45,
                      prior_halflife_var = 100,
                      prior_error_alpha = 1e-04,
                      prior_error_beta = 1e-04,
                      prior_location_gamma = NULL,
                      prior_location_range = NULL,
                      prior_max_sd_mass = 100,
                      prior_max_sd_width = 150,
                      prior_mean_pulse_count = 12,
                      sv_mass_mean = 1.5,
                      sv_width_mean = 5,
                      sv_baseline_mean = 2.6,
                      sv_halflife_mean = 45,
                      sv_error_var = 0.005,
                      sv_mass_sd = 1.6,
                      sv_width_sd = 1)


#---------------------------------------
# Fit model
#---------------------------------------
datap24<-read.delim("/Users/nichole/Repositories/FSHAnalysis2018/Data/P24.L.dat",sep=" ",header=T)
fit.24 <- fit_pulse(datap24, time = "Time", conc = "Conc",
                  spec = myspec24, iters = 50000, thin = 1, burnin = 0,
                  verbose = TRUE)


#---------------------------------------
# Check chains and posteriors
#---------------------------------------

# traceplots
bp_trace(fit.24)

# common parms posteriors
bp_posteriors(fit.24)

# pulse location posterior
bp_location_posterior(fit.24)

# Inspect pulse-specific parms
#fit$pulse_chain %>%
#  filter(iteration > 50 & pulse_num == 0) %>%
#  ggplot(aes(x = mass)) + geom_histogram()
#
#fit$pulse_chain %>%
#  filter(iteration > 50 & pulse_num == 0) %>%
#  ggplot(aes(x = iteration, y = mass)) + geom_path()


