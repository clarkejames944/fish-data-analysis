source("analysis/setup.R")

#############################################################

# function to obtain the loo results ----

comp_loo <- function(fit) {
  # extract the log likelihood from the model
  log.lik <- extract_log_lik(fit, merge_chains=FALSE)
  # obtain the relative number of effective samples
  r_eff <- relative_eff(exp(log.lik))
  # obtain the loo results
  loo(log.lik, r_eff=r_eff, cores = 4)
}

# function to obtain waic results ---
comp_waic <- function(fit) {
  # extract the log likelihood from the model
  log.lik <- extract_log_lik(fit, merge_chains=FALSE)
  # obtain the relative number of effective samples
  r_eff <- relative_eff(exp(log.lik))
  # obtain the waic results
  waic(log.lik, r_eff=r_eff, cores = 4)
}

##############################################################
## Prepare for model comparison: sex difference model ----

mod_name <- "hinge_sex_diff"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

loo_1 <- comp_loo(stan_fit)
waic_1 <- comp_waic(stan_fit)

print(loo_1)
print(waic_1)


#############################################################
## Prepare for model comparison: age difference model ----

mod_name <- "hinge_sex_age_diff"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

loo_2 <- comp_loo(stan_fit)
waic_2 <- comp_waic(stan_fit)

print(loo_2)
print(waic_2)


############################################################
## Prepare for model comparison: temperature difference model ----

mod_name <- "hinge_temp_effect"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

loo_3 <- comp_loo(stan_fit)
waic_3 <- comp_waic(stan_fit)

print(loo_3)
print(waic_3)

############################################################
## Prepare for model comparison: zone difference model ----

mod_name <- "hinge_zone_effect"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

loo_4 <- comp_loo(stan_fit)
waic_4 <- comp_waic(stan_fit)

print(loo_4)
print(waic_4)

############################################################
## Prepare for model comparison: zone and temperature difference model ----

mod_name <- "hinge_zone_temp_sex_diff"

f_loc <- sub("XX", replacement = mod_name, x = "models/XX.rds")
stan_fit <- readRDS(file = f_loc)

loo_5 <- comp_loo(stan_fit)
waic_5 <- comp_waic(stan_fit)

print(loo_5)
print(waic_5)

############################################################
## compare the models ---
#loo comparison
#1=sex_diff
#2=age effect
#3=temp effect
#4=zone effect
#5=zone and temp effect
compare(loo_1, loo_2, loo_3, loo_4, loo_5)

#waic comparison
compare(waic_1, waic_2, waic_3, waic_4, waic_5)
