source("analysis/setup.R")

#################################################################

## Read in the model: ----

# read in the model we'll use
stan_model <- readRDS(file = "models/hinge_zone_temp_sex_diff.rds")

sim_z1 <- rstan::extract(stan_model, 'sim_z1') # GLOBAL


