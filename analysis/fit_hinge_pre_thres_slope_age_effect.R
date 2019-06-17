source("analysis/setup.R")

mod_name <- "hinge_pre_thres_slope_age_effect"

###############################################################################
## Prep data for the mcmc ----

# format data for stan
fishdat_cut <- prep_stan_data(fishdat_cut)

# setting up for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###############################################################################
## Set up the initial values list ----

init_vals <- list()
# univariate parameters
init_vals$alpha        <-  0.77
init_vals$d_alpha_male <- -0.08
init_vals$d_beta_1_age <-  0.00
init_vals$beta_1       <-  0.12
init_vals$beta_2       <- -0.03
init_vals$eta          <-  0.49
init_vals$d_eta_male   <- -0.06
init_vals$sigma_fish   <-  0.05
init_vals$sigma_year   <-  0.03
init_vals$sigma        <-  0.06
# replicate the list of initial conditions to match number of chains
init_vals <- replicate(4, init_vals, simplify = FALSE)

###############################################################################
## Run the mcmc: hinge pre-threshold slope age effect model ----

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
    control = list(max_treedepth = 10, adapt_delta = 0.8)
  )
)

###############################################################################
## Save the model ----

stan_fit@stanmodel@dso <- new("cxxdso") # remove reference to stan model binary 
f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
saveRDS(stan_fit, file = f_loc)