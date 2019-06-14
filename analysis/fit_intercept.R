source("analysis/setup.R")

mod_name <- "intercept"

###############################################################################
## Prep data for the mcmc ----

# format data for stan
fishdat_cut <- prep_stan_data(fishdat_cut)

# setting up for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###############################################################################
## Run the mcmc: intercept only variation ----

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
    control = list(max_treedepth = 10, adapt_delta = 0.8)
  )
)

###############################################################################
## Save the model ----

stan_fit@stanmodel@dso <- new("cxxdso") # remove reference to stan model binary 
f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
saveRDS(stan_fit, file = f_loc)
