###############
####Time for some Stan####

#Setting up
#Libraries required
library(dplyr)
library(ggplot2)
library(rstan)
library(ggmcmc)
library(bayesplot)
library(scales)
library(coda)


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

#####The models####
#1. Simple no threshold model. -no threshold and a fixed intercept.
#2. Simple no threshold varying intercept model. -no threshold but with varying intercepts via fish individual
#3. Unhinged fixed threshold + intercept. -The discontinuous piecewise linear regression
#4. Unhinged fixed threshold, varying intercept.
#5. Unhinged fixed intercept varying threshold.
#6. Unhinged varying threshold and intercept.
#7. Hinged fixed threshold and intercept. -The continuous piecewise linear regression
#8. Hinged varying intercept fixed threshold.
#9. Hinged varying threshold fixed intercept.
#10. Hinged varying threshold and intercept.

##Running the model. 1. Simple no threshold model####
rt <- stanc(file="The Stan models/Simple no threshold model.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(fit, pars=c("intercept","beta","error","lp__", "sigma_int", "yhat[110]"))


##The extracted plot
fitted_curves <- rstan::extract(fit)
fitted_curves <- as_data_frame(fitted_curves)
head(fitted_curves)

error_hat <- median(fitted_curves$error)
yhat_hat <- rep(NA, fishdat$N_EBSm)
beta_hat <- median(fitted_curves$beta)
intercept_hat <- median(fitted_curves$intercept)

yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  yhat_hat[i] <- median(fitted_curves$yhat[,i])
}
df_post <- data.frame(list(prev=fishdat$prev,
                           oto_size=fishdat$oto_size,
                           yhat_hat=yhat_hat))

glimpse(df_post)

ggplot(df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = yhat_hat - 1.96 * error_hat,
                  ymax = yhat_hat + 1.96 * error_hat),
              fill = "lightyellow") + 
  geom_abline(aes(y=yhat_hat, slope=beta_hat, intercept=intercept_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.1)



##Plot 8 examples of the predictive models compared with the original data
simulated_data <- rstan::extract(fit)$sim_oto_size
simulated_data <- simulated_data[sample(1:4000, 8), ]

simulated_yhat <- rstan::extract(fit)$yhat
simulated_yhat <- simulated_yhat[sample(1:4000, 8), ] 

par(mfrow = c(3, 3))
rstan::plot(fishdat$prev, fishdat$oto_size)
for (i in 1:8){
  rstan::plot(fishdat$prev, simulated_data[i, ])
}

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
##Once I re-did the model it is much better

ratios_fit <- neff_ratio(fit)
print(ratios_fit)

mcmc_neff(ratios_fit, size = 2)
##Again with the number of effective sample sizes there are some that are under 0.1 (which isn't good)- 
##suggests autocorrelation for these parameters (mostly the effective sample sizes are above 0.5 though)
##Again after I re-did the model it was fine

###########################################
##Running the model. 2. Simple no threshold varying intercept model####

rt <- stanc(file="The Stan models/Simple no threshold varying intercept model.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(two_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=12, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)

pairs(two_fit, pars=c("intercept[110]","beta","error","lp__", "sigma_int", "yhat[110]"))
post <- As.mcmc.list(two_fit, pars = "yhat[110]")
plot(post)

plot(fit, pars=c("beta", "error", "yhat[110]"))

##The extracted plot
two_fitted_curves <- rstan::extract(two_fit)
two_fitted_curves <- as_data_frame(two_fitted_curves)
head(two_fitted_curves)
glimpse(two_fitted_curves)

two_error_hat <- median(two_fitted_curves$error)
two_beta_hat <- median(two_fitted_curves$beta)

two_intercept_hat <- rep(NA, fishdat$Ngroups)
for (j in 1:fishdat$Ngroups) {
  two_intercept_hat[j] <- median(two_fitted_curves$intercept[,j])
}

two_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  two_yhat_hat[i] <- median(two_fitted_curves$yhat[,i])
}

two_df_post <- data.frame(list(prev=fishdat$prev,
                           oto_size=fishdat$oto_size,
                           two_yhat_hat=two_yhat_hat))

glimpse(two_df_post)

ggplot(two_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = two_yhat_hat - 1.96 * two_error_hat,
                  ymax = two_yhat_hat + 1.96 * two_error_hat),
              fill = "lightyellow") + 
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[1]), colour = "darkred") +
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[2]), colour = "darkred") +
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[3]), colour = "darkred") +
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[4]), colour = "darkred") +
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[5]), colour = "darkred") +
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[6]), colour = "darkred") +
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[7]), colour = "darkred") +
  geom_abline(aes(y=two_yhat_hat, slope=two_beta_hat, intercept=two_intercept_hat[8]), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.1)

##Check the output
print(vs_fit)
##The 95% credible intervals seem very large for each breakpoint- not too bad now that I've constrained 
##the priors a bit but some are still fairly big

check_hmc_diagnostics(vs_fit)


#Have a look at some traceplots
posterior_vs_fit = as.array(vs_fit)

lp_vs_fit = log_posterior(vs_fit)

np_vs_fit = nuts_params(vs_fit)

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_vs_fit, pars = "bp[271]", np = np_vs_fit)
##Seems like a decent furry caterpilar plot for this parameter anyway

#Have a closer look now maybe
mcmc_trace(posterior_vs_fit, pars = "bp[271]", np = np_vs_fit, window = c(200, 350))
##Yeah looks fine really

color_scheme_set("red")
mcmc_nuts_divergence(np_vs_fit, lp_vs_fit)
##No divergence anyway

mcmc_nuts_energy(np_vs_fit)
##I think this energy plot is fine

vs_rhats = rhat(vs_fit)
print(vs_rhats)

color_scheme_set("brightblue")
mcmc_rhat(vs_rhats)
##some rhat values are above 1.1 (some seem to be around 2.67-which is really pretty bad)
##Once I re-did the model it is much better

ratios_vs_fit <- neff_ratio(vs_fit)
print(ratios_vs_fit)

mcmc_neff(ratios_vs_fit, size = 2)

###########################################################################
##Running the model. 3. Unhinged fixed threshold + intercept####

rt <- stanc(file="The Stan models/Unhinged fixed threshold + intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(three_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=17, adapt_delta=0.85)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(three_fit, pars=c("intercept1", "intercept2", "bp", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))

##The extracted plots

##The extracted plot
three_fitted_curves <- rstan::extract(three_fit)
three_fitted_curves <- as_data_frame(three_fitted_curves)
head(three_fitted_curves)
glimpse(three_fitted_curves)

three_error_hat <- median(three_fitted_curves$error)
three_beta1_hat <- median(three_fitted_curves$beta[1])
three_beta2_hat <- median(three_fitted_curves$beta[2])
three_intercept1_hat <- median(three_fitted_curves$intercept1)
three_intercept2_hat <- median(three_fitted_curves$intercept2)
three_bp_hat <- median(three_fitted_curves$bp)

three_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  three_yhat_hat[i] <- median(three_fitted_curves$yhat[,i])
}

three_df_post <- data.frame(list(prev=fishdat$prev,
                               oto_size=fishdat$oto_size,
                               three_yhat_hat=three_yhat_hat))

glimpse(three_df_post)

ggplot(three_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = three_yhat_hat - 1.96 * three_error_hat,
                  ymax = three_yhat_hat + 1.96 * three_error_hat),
              fill = "lightyellow") + 
  geom_segment(aes(x = min(fishdat$prev) , xend = three_bp_hat
                     , y = three_intercept1_hat + three_beta1_hat*min(fishdat$prev), yend = three_intercept1_hat + three_beta1_hat*three_bp_hat))+
  geom_segment(aes(x = three_bp_hat , xend = max(fishdat$prev)
                   , y = three_intercept2_hat + three_beta2_hat*three_bp_hat, yend = three_intercept2_hat + three_beta2_hat*max(fishdat$prev)))+
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.1)
 


##Check the output
print(cs_fit)
##The 95% credible intervals seem very large for each breakpoint- not too bad now that I've constrained 
##the priors a bit but some are still fairly big

check_hmc_diagnostics(cs_fit)


#Have a look at some traceplots
posterior_cs_fit = as.array(cs_fit)

lp_cs_fit = log_posterior(cs_fit)

np_cs_fit = nuts_params(cs_fit)

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cs_fit, pars = "bp[271]", np = np_cs_fit)
##Seems like a decent furry caterpilar plot for this parameter anyway

#Have a closer look now maybe
mcmc_trace(posterior_cs_fit, pars = "bp[271]", np = np_cs_fit, window = c(200, 350))
##Yeah looks fine really

color_scheme_set("red")
mcmc_nuts_divergence(np_cs_fit, lp_cs_fit)
##No divergence anyway

mcmc_nuts_energy(np_cs_fit)
##I think this energy plot is fine

cs_rhats = rhat(cs_fit)
print(cs_rhats)

color_scheme_set("brightblue")
mcmc_rhat(cs_rhats)
##some rhat values are above 1.1 (some seem to be around 2.67-which is really pretty bad)
##Once I re-did the model it is much better

ratios_cs_fit <- neff_ratio(cs_fit)
print(ratios_cs_fit)

mcmc_neff(ratios_cs_fit, size = 2)

#######################################################################
##Running the model. 4. Unhinged fixed threshold, varying intercept####

rt <- stanc(file="The Stan models/Unhinged fixed threshold, varying intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(four_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.875)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(four_fit, pars=c("intercept1[110]", "intercept2[110]", "bp", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))

##The extracted plots

four_fitted_curves <- rstan::extract(four_fit)
four_error_hat <- median(four_fitted_curves$error)
four_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  four_yhat_hat[i] <- median(four_fitted_curves$yhat[,i])
}
four_df_post <- data.frame(list(prev=fishdat$prev,
                                 oto_size=fishdat$oto_size,
                                 four_yhat_hat=four_yhat_hat))

ggplot(four_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = four_yhat_hat - 1.96 * four_error_hat,
                  ymax = four_yhat_hat + 1.96 * four_error_hat),
              fill = "lightyellow") + 
  geom_line(aes(y=four_yhat_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.05)

#######################################################################
##Running the model. 5. Unhinged fixed intercept varying threshold####

rt <- stanc(file="The Stan models/Unhinged fixed intercept varying threshold.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(five_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=15)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(five_fit, pars=c("intercept1", "intercept2", "bp[110]", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))

##The extracted plots

five_fitted_curves <- rstan::extract(five_fit)
five_error_hat <- median(five_fitted_curves$error)
five_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  five_yhat_hat[i] <- median(five_fitted_curves$yhat[,i])
}
five_df_post <- data.frame(list(prev=fishdat$prev,
                                 oto_size=fishdat$oto_size,
                                 five_yhat_hat=five_yhat_hat))

ggplot(five_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = five_yhat_hat - 1.96 * five_error_hat,
                  ymax = five_yhat_hat + 1.96 * five_error_hat),
              fill = "lightyellow") + 
  geom_line(aes(y=five_yhat_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.05)
#######################################################################
##Running the model. 6. Unhinged varying threshold and intercept####

rt <- stanc(file="The Stan models/Unhinged varying threshold and intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(six_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(six_fit, pars=c("intercept1[110]", "intercept2[110]", "bp[110]", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))

##The extracted data plot

six_fitted_curves <- rstan::extract(six_fit)
six_error_hat <- median(six_fitted_curves$error)
six_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  six_yhat_hat[i] <- median(six_fitted_curves$yhat[,i])
}
six_df_post <- data.frame(list(prev=fishdat$prev,
                                 oto_size=fishdat$oto_size,
                                 six_yhat_hat=six_yhat_hat))

ggplot(six_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = six_yhat_hat - 1.96 * six_error_hat,
                  ymax = six_yhat_hat + 1.96 * six_error_hat),
              fill = "lightyellow") + 
  geom_line(aes(y=six_yhat_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.05)
#######################################################################
##Running the model. 7. Hinged fixed threshold and intercept####

rt <- stanc(file="The Stan models/Hinged fixed threshold and intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(seven_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.825)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(seven_fit, pars=c("intercept", "bp", "mu_bp", "sigma_bp","beta1", "beta2", "error","lp__", "sigma_int", "yhat[110]", "slope_after", "intercept_after"))

##The extracted data plot

seven_fitted_curves <- rstan::extract(seven_fit)
seven_error_hat <- median(seven_fitted_curves$error)
seven_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  seven_yhat_hat[i] <- median(seven_fitted_curves$yhat[,i])
}
seven_df_post <- data.frame(list(prev=fishdat$prev,
                                 oto_size=fishdat$oto_size,
                                 seven_yhat_hat=seven_yhat_hat))

ggplot(seven_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = seven_yhat_hat - 1.96 * seven_error_hat,
                  ymax = seven_yhat_hat + 1.96 * seven_error_hat),
              fill = "lightyellow") + 
  geom_line(aes(y=seven_yhat_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.05)

#######################################################################
##Running the model. 8. Hinged varying intercept fixed threshold####

rt <- stanc(file="The Stan models/Hinged varying intercept fixed threshold.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(eight_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(eight_fit, pars=c("intercept[110]", "bp", "mu_bp", "sigma_bp","beta1", "beta2", "error","lp__", "sigma_int", "yhat[110]", "slope_after", "intercept_after[110]"))

##The extracted plots

eight_fitted_curves <- rstan::extract(eight_fit)
eight_error_hat <- median(eight_fitted_curves$error)
eight_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  eight_yhat_hat[i] <- median(eight_fitted_curves$yhat[,i])
}
eight_df_post <- data.frame(list(prev=fishdat$prev,
                                 oto_size=fishdat$oto_size,
                                 eight_yhat_hat=eight_yhat_hat))

ggplot(eight_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = eight_yhat_hat - 1.96 * eight_error_hat,
                  ymax = eight_yhat_hat + 1.96 * eight_error_hat),
              fill = "lightyellow") + 
  geom_line(aes(y=eight_yhat_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.05)
#######################################################################
##Running the model. 9. Hinged varying threshold fixed intercept####

rt <- stanc(file="The Stan models/Hinged varying threshold fixed intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(nine_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(nine_fit, pars=c("intercept", "bp[110]", "mu_bp", "sigma_bp","beta1", "beta2", "error","lp__", "sigma_int", "yhat[110]","slope_after", "intercept_after[110]"))

##The extracted plots

nine_fitted_curves <- rstan::extract(nine_fit)
nine_error_hat <- median(nine_fitted_curves$error)
nine_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  nine_yhat_hat[i] <- median(nine_fitted_curves$yhat[,i])
}
nine_df_post <- data.frame(list(prev=fishdat$prev,
                                 oto_size=fishdat$oto_size,
                                 nine_yhat_hat=nine_yhat_hat))

ggplot(nine_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = nine_yhat_hat - 1.96 * nine_error_hat,
                  ymax = nine_yhat_hat + 1.96 * nine_error_hat),
              fill = "lightyellow") + 
  geom_line(aes(y=nine_yhat_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.05)

#######################################################################
##Running the model. 10. Hinged varying threshold and intercept####

rt <- stanc(file="The Stan models/Hinged varying threshold and intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(ten_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=12)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(ten_fit, pars=c("intercept[110]", "bp[110]", "mu_bp", "sigma_bp","beta1", "beta2", "error","lp__", "sigma_int", "yhat[110]", "slope_after", "intercept_after[110]"))

##The extracted plots
ten_fitted_curves <- rstan::extract(ten_fit)
ten_error_hat <- median(ten_fitted_curves$error)
ten_yhat_hat <- rep(NA, fishdat$N_EBSm)
for (i in 1:fishdat$N_EBSm) {
  ten_yhat_hat[i] <- median(ten_fitted_curves$yhat[,i])
}
ten_df_post <- data.frame(list(prev=fishdat$prev,
                                 oto_size=fishdat$oto_size,
                                 ten_yhat_hat=ten_yhat_hat))

ggplot(ten_df_post, aes(x=prev)) +
  geom_ribbon(aes(ymin = ten_yhat_hat - 1.96 * ten_error_hat,
                  ymax = ten_yhat_hat + 1.96 * ten_error_hat),
              fill = "lightyellow") + 
  geom_line(aes(y=ten_yhat_hat), colour = "darkred") +
  geom_point(aes(y=oto_size), colour = "darkblue", alpha=0.05)
