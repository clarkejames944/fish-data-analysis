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
##Plot 8 examples of the predictive models compared with the original data
simulated_data <- rstan::extract(fit)$sim_oto_size
simulated_data <- simulated_data[sample(1:4000, 8), ]

simulated_yhat <- rstan::extract(fit)$yhat
simulated_yhat <- simulated_yhat[sample(1:4000, 8), ] 
  
ggplot(fishdat, aes(x=prev, y=oto_size))+
  geom_point()+
  geom_line(sim_yhats, aes(y=yhat))

par(mfrow = c(3, 3))
rstan::plot(fishdat$prev, fishdat$oto_size)

for (i in 1:8){
  rstan::plot(fishdat$prev, simulated_data[i, ])
}

#Make some predictive plots for general trend over original data
plot(fishdat$prev, fishdat$oto_size)


sim_yhats <- rstan::extract(fit)$yhat
yhats <- as.data.frame(sim_yhats)


df_sim_yhats <- data.frame(
  prev = fishdat$prev,
  meanGJT = apply(sim_yhats, 2, mean),
  lo80GJT = apply(sim_yhats, 2, quantile, 0.10),
  hi80GJT = apply(sim_yhats, 2, quantile, 0.90),
  lo95GJT = apply(sim_yhats, 2, quantile, 0.025),
  hi95GJT = apply(sim_yhats, 2, quantile, 0.975)
)
head(df_sim_yhats)

ggplot(df_sim_yhats,
       aes(x = prev,
           y = meanGJT)) +
  geom_ribbon(aes(ymin = lo95GJT,
                  ymax = hi95GJT),
              fill = "lightgrey") +
  geom_ribbon(aes(ymin = lo80GJT,
                  ymax = hi80GJT),
              fill = "darkgrey") +
  geom_line()

mcmc = as.matrix(fit)

## Calculate the fitted values
newdata = data.frame(x = seq(min(fishdat$prev, na.rm = TRUE), max(fishdat$prev, na.rm = TRUE),
                             len = 1000))
Xmat = model.matrix(~fishdat$prev, newdata)
coefs = mcmc[, c("beta0", "beta[1]")]
fit = coefs %*% t(Xmat)
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "HPDinterval"))


post <- As.mcmc.list(fit, pars = "oto_size[110]")
plot(post)

plot(fit, pars=c("beta", "error", "yhat[110]"))

new_y <- rstan::extract(fit, pars="yhat[1]")

ggplot(EBSm, aes(x = prev, y = oto_size)) +
  geom_point()+
  geom_line(aes(y = yhats[1]), alpha = .1) 
  scale_color_brewer(palette = "Dark2")


#check diagnostics
check_hmc_diagnostics(fit)

##Extracting and plotting
post <- rstan::extract(fit)

EBSm %>%
  data_grid(prev = seq_range(prev, n = 101)) %>%
  add_fitted_draws(post, n = 100) %>%
  ggplot(aes(x = prev, y = oto_size)) +
  geom_line(aes(y = .value), alpha = .1) +
  geom_point(data = EBSm) +
  scale_color_brewer(palette = "Dark2")

plot(fishdat$prev, fishdat$oto_size, 
     xlab = "s", ylab = "s'")
for (i in seq_along(post$lp__)) {
  segments(x0 = min(fishdat$prev), x1 = max(fishdat$prev), 
           y0 = post$intercept + post$beta * min(fishdat$prev), 
           y1 = post$intercept + post$beta * max(fishdat$prev), 
           col = alpha(3, .1))
}

ggplot(EBSm, aes(x=prev, y=oto_size))+
        geom_point()+
         xlab("s") + ylab("s'")+
        geom_line(data = post)

for (i in seq_along(post$lp__)) {
  segments(x0 = min(fishdat$prev), x1 = post$bp[i], 
           y0 = post$intercept[i] + post$beta[i, 1] * min(fishdat$prev), 
           y1 = post$intercept[i] + post$beta[i, 1] * post$bp[i], 
           col = alpha(3, .1))
  segments(x1 = max(fishdat$prev), x0 = post$bp[i], 
           y1 = post$intercept[i] + 
             post$beta[i, 2] * (max(fishdat$prev) - post$bp[i]) + 
             post$beta[i, 1] * max(fishdat$prev), 
           y0 = post$intercept[i] + post$beta[i, 1] * post$bp[i], 
           col = alpha(2, .1))
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

#fit <- stan(file="The Stan models/Simple no threshold model.stan", data = fishdat, pars = c("beta", "intercept", "error", "sigma_int","yhat"))
##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)

pairs(two_fit, pars=c("intercept[110]","beta","error","lp__", "sigma_int", "yhat[110]"))
post <- As.mcmc.list(two_fit, pars = "yhat[110]")
plot(post)

plot(fit, pars=c("beta", "error", "yhat[110]"))

new_y <- rstan::extract(fit, pars="yhat")



plot(fishdat$prev, fishdat$oto_size, 
     xlab = 's', ylab = 's following')
for (i in seq_along(vs_post$lp__)) {
  segments(x0 = min(fishdat$prev), x1 = vs_post$bp[i], 
           y0 = vs_post$intercept[i] + vs_post$beta[i, 1] * min(fishdat$prev), 
           y1 = vs_post$intercept[i] + vs_post$beta[i, 1] * vs_post$bp[i], 
           col = alpha(3, .1))
  segments(x1 = max(fishdat$prev), x0 = vs_post$bp[i], 
           y1 = vs_post$intercept[i] + 
             vs_post$beta[i, 2] * (max(fishdat$prev) - vs_post$bp[i]) + 
             vs_post$beta[i, 1] * max(fishdat$prev), 
           y0 = vs_post$intercept[i] + vs_post$beta[i, 1] * vs_post$bp[i], 
           col = alpha(2, .1))
}



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
system.time(three_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=12, adapt_delta=0.85)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(three_fit, pars=c("intercept1", "intercept2", "bp", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))
print(three_fit, pars=c("intercept1", "intercept2", "bp", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))
three_post <- rstan::extract(three_fit) 

rstan::plot(fishdat$prev, fishdat$oto_size, 
     xlab = 's', ylab = 's following')
 segments(x0 = min(fishdat$prev), x1 = three_post$bp, 
           y0 = three_post$intercept1 + three_post$beta[1] * min(fishdat$prev), 
           y1 = three_post$intercept1 + three_post$beta[1] * three_post$bp,
          col=4)
 segments(x1 = max(fishdat$prev), x0 = three_post$bp, 
           y1 = three_post$intercept2 + three_post$beta[2] * (max(fishdat$prev)), 
           y0 = three_post$intercept2 + three_post$beta[2] * three_post$bp,
          col = 5)


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


#######################################################################
##Running the model. 5. Unhinged fixed intercept varying threshold####

rt <- stanc(file="The Stan models/Unhinged fixed intercept varying threshold.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(five_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=15)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(five_fit, pars=c("intercept1", "intercept2", "bp[110]", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))


#######################################################################
##Running the model. 6. Unhinged varying threshold and intercept####

rt <- stanc(file="The Stan models/Unhinged varying threshold and intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(six_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(six_fit, pars=c("intercept1[110]", "intercept2[110]", "bp[110]", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int1", "sigma_int2", "yhat[110]"))


#######################################################################
##Running the model. 7. Hinged fixed threshold and intercept####

rt <- stanc(file="The Stan models/Hinged fixed threshold and intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(seven_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(seven_fit, pars=c("intercept1", "intercept2", "bp", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int", "yhat[110]"))


#######################################################################
##Running the model. 8. Hinged varying intercept fixed threshold####

rt <- stanc(file="The Stan models/Hinged varying intercept fixed threshold.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(eight_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(eight_fit, pars=c("intercept1[110]", "intercept2[110]", "bp", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int", "yhat[110]"))


#######################################################################
##Running the model. 9. Hinged varying threshold fixed intercept####

rt <- stanc(file="The Stan models/Hinged varying threshold fixed intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(nine_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(nine_fit, pars=c("intercept1", "intercept2", "bp[110]", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int", "yhat[110]"))


#######################################################################
##Running the model. 10. Hinged varying threshold and intercept####

rt <- stanc(file="The Stan models/Hinged varying threshold and intercept.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(ten_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(ten_fit, pars=c("intercept1[110]", "intercept2[110]", "bp[110]", "mu_bp", "sigma_bp","beta[1]", "beta[2]", "error","lp__", "sigma_int", "yhat[110]"))

