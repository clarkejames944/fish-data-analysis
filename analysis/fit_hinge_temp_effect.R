source("analysis/setup.R")

mod_name <- "hinge_temp_effect"

###############################################################################
## Prep data for the mcmc ----

# format data for stan
fishdat_cut <- prep_stan_data(fishdat_cut)

# setting up for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###############################################################################
## Set up the initial values list ----

# read in the model we'll take initial values from
init_model <- readRDS(file = "models/hinge_sex_diff.rds")
# calculate posterior mean of univariate parameters
pars <- 
  c("alpha", "d_alpha_male", "d_eta_male",
    "d_beta_1_male", "d_beta_2_male",
    "beta_1", "beta_2", "eta", "sigma_fish", "sigma_year", "sigma")
init_vals <- lapply(rstan::extract(init_model, pars), mean)
# calculate posterior mean of vector-valued parameters
init_vals$u_fish <- apply(rstan::extract(init_model, "u_fish")[[1]], 2, mean)
init_vals$u_year <- apply(rstan::extract(init_model, "u_year")[[1]], 2, mean)
# add any additional parameters in new model
init_vals$d_alpha_temp <- 0.0 
init_vals$d_eta_temp <- 0.0 
init_vals$d_beta_1_temp <- 0.0 
init_vals$d_beta_2_temp <- 0.0 
# replicate the list of initial conditions to match number of chains
init_vals <- replicate(4, init_vals, simplify = FALSE)
# remove the massive object
rm(init_model); gc()


###############################################################################
## Run the mcmc: hinge temperature effect model ----

# compile model
f_loc <- sub("XX", replacement = mod_name, x = "stan/XX.stan")
rt <- stanc(file = f_loc)
sm <- stan_model(stanc_ret = rt, verbose = FALSE)

# 
system.time(
  stan_fit <- sampling(
    sm,
    data = fishdat_cut,
    seed = 1,
    iter = 2000,
    chains = 4,
    init = init_vals,
    control = list(max_treedepth = 11, adapt_delta = 0.8)
  )
)

###############################################################################
## Save the model ----

stan_fit@stanmodel@dso <- new("cxxdso") # remove reference to stan model binary 
f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
saveRDS(stan_fit, file = f_loc)