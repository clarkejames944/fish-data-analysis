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

###############################################################################
## Summarise the model: hinge sex diff ----

mod_name <- "hinge_sex_diff"

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
    "d_beta_1_male", "d_beta_2_male",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge intercept age effect ----

mod_name <- "hinge_intercept_age_effect"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_alpha_age",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge threshold age effect ----

mod_name <- "hinge_threshold_age_effect"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_eta_age",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge pre-threshold slope age effect ----

mod_name <- "hinge_pre_thres_slope_age_effect"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_beta_1_age",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()
###############################################################################
## Summarise the model: hinge intercept and threshold age effect ----

mod_name <- "hinge_intercept+threshold_age_effect"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_alpha_age", "d_eta_age",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge intercept and pre-threshold slope age effect ----

mod_name <- "hinge_intercept+pre_thres_slope_age_effect"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_alpha_age", "d_beta_1_age",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge threshold and pre-threshold slope age effect ----

mod_name <- "hinge_threshold+pre_thres_slope_age_effect"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_eta_age", "d_beta_1_age",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge all fixed effects age effect ----

mod_name <- "hinge_age_effect_all_fixed_effects"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_alpha_age",
    "d_eta_age", "d_beta_1_age",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge sex and age effect ----

mod_name <- "hinge_sex_age_diff"

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
    "d_beta_1_male", "d_beta_2_male",
    "d_alpha_age",
    "d_eta_age", "d_beta_1_age",
    "d_beta_2_age",
    "eta", "d_eta_male",
    "d_eta_age",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge temperature effect ----

mod_name <- "hinge_temp_effect"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

# pairs plots
f_loc <- sub("XX", replacement = mod_name, x = "figures/XX_pairs.pdf")
pdf(file = f_loc, w = 15, h = 15)
pairs(
  stan_fit,
  pars = c(
    "alpha", "d_alpha_male", "d_alpha_temp",
    "beta_1", "beta_2", 
    "d_beta_1_male", "d_beta_2_male", "d_beta_1_temp",
    "d_beta_2_temp",
    "d_eta_temp",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()

###############################################################################
## Summarise the model: hinge zone effect ----

mod_name <- "hinge_zone_effect"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

# pairs plots
f_loc <- sub("XX", replacement = mod_name, x = "figures/XX_pairs.pdf")
pdf(file = f_loc, w = 15, h = 15)
pairs(
  stan_fit,
  pars = c(
    "alpha", "d_alpha_male", "d_alpha_zone",
    "beta_1", "beta_2", 
    "d_beta_1_male", "d_beta_2_male", "d_beta_1_zone",
    "d_beta_2_zone",
    "d_eta_zone",
    "eta", "d_eta_male",
    "sigma_fish", "sigma_year", "sigma",
    "lp__"
  )
)
dev.off()
