###############
####Time for some Stan####

#Setting up
#Libraries required
library(dplyr)
library(ggplot2)
library(rstan)
library(rstanarm)
library(lattice)
library(reshape2)
library(ggmcmc)
library(bayesplot)

#Data required
fish <- read.csv("https://raw.githubusercontent.com/clarkejames944/fish-data-analysis/master/otoliths%20(working)/data_derived/data_otolith_complete.csv")

##Remove fish with unknown gender from the dataset
fish <- fish %>% filter(sex !='U')

##Must create otoloth size information. Federico did this by;
oto_size <- rep(0,23625)
for (i in 1:23625) oto_size[i] <- sum(fish$Increment[((i-fish$Age[i])+1):i])
fish<-data.frame(cbind(fish,oto_size))

##Then create log of otolith size
fish<-mutate(fish, log_oto_size=log(oto_size))
write.csv(fish,"extended_data_oto.csv")

##Creation of a GAM of initial otolith size against otolith size in 2nd year
##Need to create a new variable of oto_size in of year 1 and not year 1

##not_1 is a new data frame containing all measures not equal to year 1 fish
##Following Fede's method of doing this
fish <- read.csv("extended_data_oto.csv")
attach(fish)
ind_not_1 <- which(Age!=1)
not_1 <- filter(fish, Age!=1 )
not_1 <- mutate(not_1, prev=fish[ind_not_1-1,33])
write.csv(not_1,"parallel_data_oto.csv")

#create subsets
EBSm <- not_1 %>% filter(sex=='M', zone=='EBS')
EBSf <- not_1 %>% filter(sex=='F', zone=='EBS')
ETASm <- not_1 %>% filter(sex=='M', zone=='ETAS')
ETASf <- not_1 %>% filter(sex=='F', zone=='ETAS')
NSWm <- not_1 %>% filter(sex=='M', zone=='NSW')
NSWf <- not_1 %>% filter(sex=='F', zone=='NSW')
WTASm <- not_1 %>% filter(sex=='M', zone=='WTAS')
WTASf <- not_1 %>% filter(sex=='F', zone=='WTAS')


##Setting up for Stan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(1)

##Specifying the data required for this model

fishdat <- list(N_EBSm=nrow(EBSm),
                Ngroups = length(unique(EBSm$FishID)),
                oto_size=(EBSm$oto_size)^2,
                Age= EBSm$Age,
                fishID= as.numeric(factor(EBSm$FishID)),
                prev=(EBSm$prev)^2
)

##Running the model
rt <- stanc(file="Stan code.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
#system.time(fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=15)))



##Check the output
print(fit)
##The 95% credible intervals seem very large for each breakpoint- not too bad now that I've constrained 
##the priors a bit but some are still fairly big


#shinystan::launch_shinystan(fit)

check_treedepth(fit)

##All the 4000 iterations the maximum tree depth of 10 but zero saturated the maximum of 15

check_energy(fit)
##The E-BFMI seems okay

check_divergences(fit)
##0 of 4000 iterations ended with a divergence

check_hmc_diagnostics(fit)


#Have a look at some traceplots
posterior_fit = as.array(fit)

lp_fit = log_posterior(fit)

np_fit = nuts_params(fit)

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_fit, pars = "bp[271]", np = np_fit)
##Seems like a decent furry caterpilar plot for this parameter anyway

#Have a closer look now maybe
mcmc_trace(posterior_fit, pars = "bp[271]", np = np_fit, window = c(200, 350))
##Yeah looks fine really

color_scheme_set("red")
mcmc_nuts_divergence(np_fit, lp_fit)
##No divergence anyway

mcmc_nuts_energy(np_fit)
##I think this energy plot is fine

rhats = rhat(fit)
print(rhats)

color_scheme_set("brightblue")
mcmc_rhat(rhats)
##some rhat values are above 1.1 (some seem to be around 2.67-which is really pretty bad)

ratios_fit <- neff_ratio(fit)
print(ratios_fit)

mcmc_neff(ratios_fit, size = 2)
##Again with the number of effective sample sizes there are some that are under 0.1 (which isn't good)- 
##suggests autocorrelation for these parameters (mostly the effective sample sizes are above 0.5 though)

color_scheme_set("red")
ppc_dens_overlay(y = fit$oto_size[1],
                 yrep = posterior_predict(fit, draws = 50))

