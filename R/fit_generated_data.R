#-------------------------------------------------------------------------------
# generate_data.R
#   Generate simulated data for debugging single subject model
#
# Author: M.Mulvahill
# Notes:
#   - generating 4 series, 3 clean and 1 noisy
# 
#   - Other than changing mass_mean to 1.25 and mass_sd to 0.5 in the noisy
#   series, the values for simulation are defaults:
#   baseline              = 2.6
#   mean pulse mass       = 3.5
#   SD/var of pulse mass  = 1.6 (sd)
#   mean pulse width      = 40
#   SD/var of pulse width = 20 (sd)
#   half-life             = 45
#   model error           = 0.005
# 
#-------------------------------------------------------------------------------

setwd("~/Projects/BayesPulse/Software/singlesubject_debug")

if (!require(devtools))  install.packages("devtools")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(furrr)) install.packages("furrr")
if (!require(pulsatile)) devtools::install_github("bayespulse/pulsatile")
library(tidyverse)
library(pulsatile)

# Load data
load("./data/sim_series_r_versions.RData")

# Fit clean series
set.seed(2018-09-06)
plan(multiprocess)

specs <- map(1:3, ~ pulse_spec(location_prior_type = "strauss",
                             prior_location_gamma = 0,
                             prior_location_range = 30))
specs <- c(specs, list(pulse_spec(location_prior_type = "strauss",
                             prior_location_gamma = 0,
                             prior_location_range = 30,
                             prior_mass_mean = 1.5, sv_mass_mean = 1.5)))

fits <- future_map2(list(clean_series_1, clean_series_2, clean_series_3,
                         one_noisy_series),
                    specs, ~ fit_pulse(.x, spec = .y))


map(fits, bp_trace)






