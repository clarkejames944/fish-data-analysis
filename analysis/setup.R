#######################################################################
## Load and attach libraries ----

library(dplyr)
library(ggplot2)
library(rstan)
library(ggmcmc)
library(bayesplot)
library(scales)
library(coda)
library(stringr)
library(mgcv)
library(loo)
library(car)
library(Matrix)
library(viridis)
library(reshape2)

#######################################################################
## Function definitions ----

# format the data for the stan model
prep_stan_data <- function(data) {
  list(
    n_obs   = nrow(data),
    n_fish  = length(unique(data$FishID)),
    n_year  = length(unique(data$Year)),
    id_fish = as.integer(factor(data$FishID)),
    id_year = as.integer(factor(data$Year)),
    a       = as.integer(data$Age),
    temp    = as.numeric(data$bottomtemp1),
    is_EBS  = as.integer(data$zone == "EBS"),
    is_ETAS = as.integer(data$zone == "ETAS"),
    is_WTAS = as.integer(data$zone == "WTAS"),
    is_NSW  = as.integer(data$zone == "NSW"),
    is_f    = as.integer(data$sex == "F"),
    is_m    = as.integer(data$sex == "M"),
    z0      = data$z0^2,
    z1      = data$z1^2
  )
}

#######################################################################
## Read in data and derive variables ----

repo <- "https://raw.githubusercontent.com/clarkejames944"
floc <- "fish-data-analysis/master/otoliths%20(working)/data_derived"
fnme <- "data_otolith_complete.csv"

# Read in the data 
fish <- read.csv(file.path(repo, floc, fnme)) 

# Create new individual-level variables
fish <- fish %>% 
  group_by(FishID) %>% 
  mutate(
    # create otoloth size variable (== cumulative increment)
    z0 = cumsum(Increment),
    # otolith size next year
    z1 = lead(z0)
   ) %>% 
  ungroup()

# Create factor versions of variables
fish <- fish %>% 
  mutate(
    # create factor version of age
    Age_f = factor(Age),
    # create factor version of year
    Year_f = factor(Year)
  )

# select useable cases
fish <- fish %>% 
  filter(
    # remove fish with unknown sex from the dataset
    sex != 'U',
    # remove the last obs per fish
    !is.na(z1)
  )

#######################################################################
## Create subset we want to work with ----

# filter the data for 
fishdat_cut <- fish %>% 
  # select obs we want to work with 
  filter(
    # trim the range of ages included
    Age >= 2, Age <= 8
  )


