source("analysis/setup.R")

#################################################################

## Read in the model: ----

# read in the model we'll use
stan_model <- readRDS(file = "models/hinge_zone_temp_sex_diff.rds")


## extract the simulated values
sim_z1 <- rstan::extract(stan_model, 'sim_z1') 
sim_z1 <- as.matrix(stan_model, pars= "sim_z1")



##### check the overlay of simulated y values against real y values
ppc_dens_overlay(fishdat_cut$z1^2, sim_z1[1:500, ])


#### compare the mean of y against mean of iterations
ppc_stat(y = fishdat_cut$z1^2, yrep = sim_z1, stat = "mean")

#### compare the mean of y against mean of iterations but grouped by age
ppc_stat_grouped(y = fishdat_cut$z1^2, yrep = sim_z1, group= fishdat_cut$Age, stat = "mean")

#### same thing as above but frequency polygon
ppc_stat_freqpoly_grouped(y = fishdat_cut$z1^2, yrep = sim_z1, group= fishdat_cut$Age, stat = "mean")
#### same thing as above but scatter plot
ppc_scatter_avg_grouped(y = fishdat_cut$z1^2, yrep = sim_z1, group= fishdat_cut$Age)



### the following don't seem to work properly- ok they do now
ppc_freqpoly_grouped(y = fishdat_cut$z1^2, yrep = sim_z1[1:8, ], group= fishdat_cut$Age)

ppc_intervals_grouped(y = fishdat_cut$z1^2, yrep = sim_z1[1:8, ], group= fishdat_cut$Age)

ppc_ribbon_grouped(y = fishdat_cut$z1^2, yrep = sim_z1[1:20, ], group= fishdat_cut$Age)

ppc_error_hist_grouped(y = fishdat_cut$z1^2, yrep = sim_z1[1:10, ], group= fishdat_cut$Age)

ggplot(fishdat_cut, aes(x=z1))+
  geom_histogram()+
  ppc_dens_overlay(fishdat_cut$z1^2, sim_z1[1:500, ])

available_ppc()
