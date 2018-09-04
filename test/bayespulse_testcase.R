################################################################################
# The pulsatile package's mean_contrib and mean_conc nails this series
################################################################################


load("../experiments_libpulsatile/lib/libpulsatile/R-package/data-raw/test_data.RData")

test_data
library(devtools)
document()
install("./")
library(pulsatile)
library(tidyverse)

myspec <- pulse_spec(location_prior_type = "strauss",
                     prior_mass_mean = 3.5,
                     prior_mass_var = 100,
                     prior_width_mean = 42,
                     prior_width_var = 100,
                     prior_baseline_mean = 2.6,
                     prior_baseline_var = 100,
                     prior_halflife_mean = 45,
                     prior_halflife_var = 100,
                     prior_error_alpha = 1e-04,
                     prior_error_beta = 1e-04,
                     prior_location_gamma = 0,
                     prior_location_range = 40,
                     prior_max_sd_mass = 10,
                     prior_max_sd_width = 150,
                     prior_mean_pulse_count = 12,
                     sv_mass_mean = 3.5,
                     sv_width_mean = 42,
                     sv_baseline_mean = 2.6,
                     sv_halflife_mean = 45,
                     sv_error_var = 0.05,
                     sv_mass_sd = 10,
                     sv_width_sd = 10,
                     pv_baseline = 0.5,
                     pv_halflife = 45,
                     pv_mean_pulse_mass = 2,
                     pv_mean_pulse_width = 5,
                     pv_indiv_pulse_mass = 2,
                     pv_indiv_pulse_width = 2,
                     pv_sd_pulse_mass = 2,
                     pv_sd_pulse_width = 10,
                     pv_sdscale_pulse_mass = 1,
                     pv_sdscale_pulse_width = 1,
                     pv_pulse_location = 10)

fit <- fit_pulse(test_data, spec = myspec, iters = 250000, thin = 50, burnin = 25000,
                 verbose = TRUE)

fit$pulse_chain %>% 
  group_by(iteration) %>%
  summarise(mean_eta_mass = mean(eta_mass),
            mean_eta_width = mean(eta_width),
            mean_mass = mean(mass),
            mean_width = mean(width),
            sd_width = sd(width)) %>%
  summarise(mean_eta_mass = mean(mean_eta_mass),
            mean_eta_width = mean(mean_eta_width), 
            mean_mass = mean(mean_mass),
            mean_width = mean(mean_width),
            mean_sd_width = mean(sd_width))


# Initial mean_concentration
pulsatile_mean_conc <- 
  c(1.013895, 1.266416, 1.607703, 1.747340, 1.715250, 1.640688, 1.567391,
    1.499826, 1.438032, 1.381839, 1.331019, 1.285298, 1.244371, 1.208004,
    1.177663, 1.165266, 1.193239, 1.240540, 1.252347, 1.228081, 1.194994,
    1.164199, 1.136945, 1.112973, 1.091955, 1.073580, 1.057555, 1.044955,
    1.131759, 1.609830, 1.803182, 1.736639, 1.656223, 1.583179, 1.576934,
    1.843972, 2.044371, 1.996193, 1.901049, 1.809710, 1.724157, 1.644498,
    1.570765, 1.502919, 1.440853, 1.384397, 1.333327, 1.288133, 1.431296,
    1.973725, 1.950407, 1.856336, 1.767763, 1.685041, 1.608238, 1.537350,
    1.472305, 1.412967, 1.359137, 1.310566, 1.266965, 1.228012, 1.193367,
    1.162682, 1.135606, 1.111799, 1.091975, 1.157628, 1.607757, 1.817146,
    1.753454, 1.671904, 1.596084, 1.526570, 1.479363, 1.555295, 1.703273,
    1.703314, 1.631954, 1.559318, 1.492472, 1.436823, 1.494401, 1.835416,
    2.043448, 2.001371, 1.906906, 1.815230, 1.729312, 1.649285, 1.575183, 
    1.506977, 1.690332, 1.613409, 1.542109, 1.476686, 1.468308, 1.809871, 
    1.776511, 1.693237, 1.615827, 1.544335, 1.478698, 1.418784, 1.364401,
    1.315304, 1.271208, 1.231795, 1.196726, 1.165651, 1.138222, 1.114094,
    1.092937, 1.074437, 1.058315, 1.044656, 1.037847, 1.064495, 1.176554,
    1.349225, 1.464914, 1.479030, 1.436857, 1.383045, 1.332254, 1.286411,
    1.245363, 1.208781, 1.176319, 1.147823, 1.140473, 1.306750, 1.526503,
    1.522298, 1.460701, 1.402427, 1.349607, 1.301994, 1.259292, 1.221176,
    1.187303, 1.157324, 1.130888, 1.107657)

# Initial likelihood
pulsatile_likelihood <- 74.693770


pulsatile_mean_conc %>%
  exp %>%
  as_data_frame %>%
  rename(mean_concentration = value) %>%
  mutate(time = 1:n()*10) %>%
  full_join(test_data$data) %>%
  ggplot(aes(x = time)) + 
  geom_path(aes(y = mean_concentration)) +
  geom_path(aes(y = concentration))


likelihood_fn <- function(mean_conc, conc) {

  conc <- log(conc)

  sigma <- 0.050000
  lsigma <- log(sigma)
  x = 0
  for (i in 1:length(conc)) {
    x <- x + (conc[i] - mean_conc[i]) * (conc[i] - mean_conc[i])
  }
  x <- x / (-2.0 * sigma)
  x <- x + (-0.5 * length(conc) * (1.8378771 + lsigma))

  return(x)

}
likelihood_fn(pulsatile_mean_conc, test_data$data$concentration)

