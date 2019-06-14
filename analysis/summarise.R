source("analysis/setup.R")

###############################################################################
## Summarise the model: only intercept variation ----

mod_name <- "intercept"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

# pairs plots
f_loc <- sub("XX", replacement = mod_name, x = "figures/XX_pairs.pdf")
pdf(file = f_loc, w = 15, h = 15)
pairs(
  stan_fit,
  pars = c(
    "alpha", "d_alpha_male", "d_alpha_mature",
    "beta_1", "beta_2",
    "eta",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: intercept + threshold variation ----

mod_name <- "intercept_and_threshold"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

# pairs plots
f_loc <- sub("XX", replacement = mod_name, x = "figures/XX_pairs.pdf")
pdf(file = f_loc, w = 15, h = 15)
pairs(
  stan_fit,
  pars = c(
    "alpha", "d_alpha_male", "d_alpha_mature",
    "beta_1", "beta_2", 
    "gamma", "eta",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge basemodel ----

mod_name <- "hinge_basemodel"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

# pairs plots
f_loc <- sub("XX", replacement = mod_name, x = "figures/XX_pairs.pdf")
pdf(file = f_loc, w = 15, h = 15)
pairs(
  stan_fit,
  pars = c(
    "alpha", "d_alpha_male",
    "beta_1", "beta_2", 
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()


