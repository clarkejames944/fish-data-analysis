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
library(stringr)
library(mgcv)
library(loo)


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
##I am following Fede's method of doing this
fish <- read.csv("extended_data_oto.csv")
attach(fish)
ind_not_1 <- which(Age!=1)
not_1 <- filter(fish, Age!=1 )
not_1 <- mutate(not_1, prev=fish[ind_not_1-1,33])
write.csv(not_1,"parallel_data_oto.csv")

#create subsets

#Male subset
M <- not_1 %>% filter(sex=='M')

#Zone and Sex subsets
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

##Specifying the data required for the stan model
##I am squaring the variables after the exploratory analysis

fishdat <- list(N_EBSm=nrow(EBSm),
                Ngroups = length(unique(EBSm$FishID)),
                oto_size=(EBSm$oto_size)^2,
                Age= EBSm$Age,
                fishID= as.numeric(factor(EBSm$FishID)),
                prev=(EBSm$prev)^2
            )




##After finding the wrong breakpoint in the model multiple times we decided to cut the dataset after this
##wrong breakpoint (just for exploratory reasons)

#Cut off the first wrong breakpoint
EBSm_cut <- filter(EBSm, prev>0.5)

##update the data for the stan model
fishdat_cut <- list(N_EBSm=nrow(EBSm_cut),
                    Ngroups = length(unique(EBSm_cut$FishID)),
                    oto_size=(EBSm_cut$oto_size)^2,
                    Age= EBSm_cut$Age,
                    fishID= as.numeric(factor(EBSm_cut$FishID)),
                    prev=(EBSm_cut$prev)^2
)



fishdat_male <- list(N_EBSm=nrow(M),
                Ngroups = length(unique(M$FishID)),
                oto_size=(M$oto_size)^2,
                Age= M$Age,
                fishID= as.numeric(factor(M$FishID)),
                prev=(M$prev)^2
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

##The varying slope unhinged models
#11. Unhinged fixed threshold + intercept (varying slopes).
#12. Unhinged fixed threshold, varying intercept (varying slopes).
#13. Unhinged fixed intercept varying threshold (varying slopes).
#14. Unhinged varying threshold and intercept (varying slopes).

##Running the model. 1. Simple no threshold model####
rt <- stanc(file="The Stan models/Simple no threshold model.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
system.time(fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

system.time(fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))
##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
#Will give an idea of whether the parameters are being sampled well
pairs(fit, pars=c("alpha","beta","epsilon","lp__", "yhat[110]"))

##Extract the summary from the fit to plot the predicted lines over the original dataset
#Set-up as a data-frame
summary1 <- summary(fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  dplyr::select(variable, everything()) %>% 
  as_data_frame()

#have a look at what we have now
head(summary1)

##select only the simulated data from this new dataset
new_oto_size <- summary1 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

##take another look at what we have in terms of the credible intervals and median of these predictions
head(new_oto_size)

#Create a new data frame in which the mean and credible intervals are present for each data point
fishdat1 <- EBSm %>% mutate(mean_oto_size = new_oto_size$mean,
                               lower = new_oto_size$`2.5%`,
                               upper = new_oto_size$`97.5%`,
                               pred_dif = ((oto_size)^2 - new_oto_size2$mean))


fishdat1a <- fishdat1 %>% filter(pred_dif>0)
fishdat1b <- fishdat1 %>% filter(pred_dif<0)

high1 <- dplyr::summarise(fishdat1a, avg= mean(pred_dif))
low1 <- dplyr::summarise(fishdat1b, avg= mean(pred_dif))
print(high1)
print(low1)

q1 <- as.vector(quantile(fishdat1$pred_dif))

##Create the plot over the original data with the simulated values
fishdat1 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour="gold2", alpha=0.25, size=0.25) +
  geom_ribbon(aes(x=(prev)^2, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_classic()


#Graph of observed against predicted values for each
fishdat1 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat1 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept = 0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik1 <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = TRUE)
loo1 <- loo(lik1, save_psis = TRUE)
waic1 <- waic(lik1)


##Check the output
#shinystan::launch_shinystan(fit)
check_hmc_diagnostics(fit)


#Have a look at some traceplots
posterior_fit = as.array(fit)

lp_fit = log_posterior(fit)

np_fit = nuts_params(fit)

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_fit, pars = "yhat[271]", np = np_fit)

mcmc_trace(posterior_fit, pars = "yhat[271]", np = np_fit, window = c(200, 350))

color_scheme_set("red")
mcmc_nuts_divergence(np_fit, lp_fit)

mcmc_nuts_energy(np_fit)

rhats = rhat(fit)
print(rhats)

color_scheme_set("brightblue")
mcmc_rhat(rhats)

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

system.time(two_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=13, adapt_delta=0.8)))

system.time(two_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=13, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)

pairs(two_fit, pars=c("alpha[110]","beta","epsilon","lp__", "sigma_alpha", "yhat[110]"))


##Extract the summary from the fit to plot the predicted lines
summary2 <- summary(two_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  dplyr::select(variable, everything()) %>% 
  as_data_frame()

#Have a glance at what the predicted data is like
head(summary2)

#Only want the sim_oto_size parameter for the predictions
new_oto_size2 <- summary2 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))


head(new_oto_size2)


fishdat2 <- EBSm %>% mutate(mean_oto_size = new_oto_size2$mean,
                            lower = new_oto_size2$`2.5%`,
                            upper = new_oto_size2$`97.5%`,
                            pred_dif = ((oto_size)^2 - new_oto_size2$mean))

fishdat2a <- fishdat2 %>% filter(pred_dif>0)
fishdat2b <- fishdat2 %>% filter(pred_dif<0)

high2 <- dplyr::summarise(fishdat2a, avg= mean(pred_dif))
low2 <- dplyr::summarise(fishdat2b, avg= mean(pred_dif))
print(high2)
print(low2)


fishdat2 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()

  geom_ribbon(aes(x= (prev)^2, ymin = lower, ymax = upper), alpha = 0.25)

#Graph of observed against predicted values for each
fishdat2 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()
  
#Graph of oto-size against difference between predicted and observed data
  fishdat2 %>% 
    ggplot() +
    geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
    geom_hline(yintercept=0)+
    theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik2 <- extract_log_lik(two_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo2 <- loo(lik2, save_psis = TRUE)
waic2 <- waic(lik2)
print(loo2)

#Compare these first two models
loo::compare(loo1,loo2)  
#Positive value of output shows that the second model is the better for describing the data

q2 <- as.vector(quantile(fishdat2$pred_dif))
            
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

#For the original dataset
system.time(three_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=15, adapt_delta=0.8)))

#For the cut dataset
system.time(three_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

#For the complete male dataset
system.time(three_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))



##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(three_fit, pars=c("alpha1", "alpha2", "beta", "epsilon","lp__", "eta"))

##Extract the summary from the fit to plot the predicted lines
summary3 <- summary(three_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary3)

new_oto_size3 <- summary3 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size3)

fishdat3 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size3$mean,
                            lower = new_oto_size3$`2.5%`,
                            upper = new_oto_size3$`97.5%`,
                            pred_dif = ((oto_size)^2 - new_oto_size3$mean))

fishdat3a <- fishdat3 %>% filter(pred_dif>0)
fishdat3b <- fishdat3 %>% filter(pred_dif<0)

high3 <- dplyr::summarise(fishdat3a, avg= mean(pred_dif))
low3 <- dplyr::summarise(fishdat3b, avg= mean(pred_dif))
print(high3)
print(low3)


fishdat3 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  geom_ribbon(aes(x=(prev)^2, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_classic()



#Graph of observed against predicted values for each
fishdat3 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat3 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik3 <- extract_log_lik(three_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo3 <- loo(lik3, save_psis = TRUE)
waic3 <- waic(lik3)
print(loo3)
plot(loo3)


q3 <- as.vector(quantile(fishdat3$pred_dif))

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

#For the main dataset
system.time(four_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

#For the cut dataset
system.time(four_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

#For the complete male dataset
system.time(four_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(four_fit, pars=c("alpha1[110]", "alpha2[110]", "eta","beta", "epsilon","lp__", "sigma_alpha", "yhat[110]"))


##Extract the summary from the fit to plot the predicted lines
summary4 <- summary(four_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary4)

new_oto_size4 <- summary4 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size4)

fishdat4 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size4$mean,
                                lower = new_oto_size4$`2.5%`,
                                upper = new_oto_size4$`97.5%`,
                                pred_dif = ((oto_size)^2 - new_oto_size4$mean))

fishdat4a <- fishdat4 %>% filter(pred_dif>0)
fishdat4b <- fishdat4 %>% filter(pred_dif<0)

high4 <- dplyr::summarise(fishdat4a, avg= mean(pred_dif))
low4 <- dplyr::summarise(fishdat4b, avg= mean(pred_dif))
print(high4)
print(low4)

q4 <- as.vector(quantile(fishdat4$pred_dif))

fishdat4 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()

#Graph of observed against predicted values for each
fishdat4 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat4 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik4 <- extract_log_lik(four_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo4 <- loo(lik4, save_psis = TRUE)
loo4 <- loo(lik4, k_threshold= 0.7)
waic4 <- waic(lik4)
print(loo4)
plot(loo4)

loo::compare(loo4,loo14)

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

#For the original dataset
system.time(five_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the cut dataset
system.time(five_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=11)))

#For the complete male dataset
system.time(five_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(five_fit, pars=c("alpha1", "alpha2", "eta[110]", "mu_eta", "sigma_eta","beta", "epsilon","lp__", "yhat[110]"))


##Extract the summary from the fit to plot the predicted lines
summary5 <- summary(five_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary5)

new_oto_size5 <- summary5 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size5)

fishdat5 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size5$mean,
                                lower = new_oto_size5$`2.5%`,
                                upper = new_oto_size5$`97.5%`,
                                pred_dif = ((oto_size)^2 - new_oto_size5$mean))

fishdat5a <- fishdat5 %>% filter(pred_dif>0)
fishdat5b <- fishdat5 %>% filter(pred_dif<0)

high5 <- dplyr::summarise(fishdat5a, avg= mean(pred_dif))
low5 <- dplyr::summarise(fishdat5b, avg= mean(pred_dif))
print(high5)
print(low5)

q5 <- as.vector(quantile(fishdat5$pred_dif))

fishdat5 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()


#Graph of observed against predicted values for each
fishdat5 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat5 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik5 <- extract_log_lik(five_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo5 <- loo(lik5, save_psis = TRUE)
waic5 <- waic(lik5)
print(loo5)
plot(loo5)

loo::compare(loo4,loo5)

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

#For original dataset
system.time(six_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the cut dataset
system.time(six_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the complete male dataset
system.time(six_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(six_fit, pars=c("alpha1[110]", "alpha2[110]", "eta[110]", "mu_eta", "sigma_eta","beta", "epsilon","lp__", "sigma_alpha", "yhat[110]"))


##Extract the summary from the fit to plot the predicted lines
summary6 <- summary(six_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary6)

new_oto_size6 <- summary6 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size6)

fishdat6 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size6$mean,
                                lower = new_oto_size6$`2.5%`,
                                upper = new_oto_size6$`97.5%`,
                                pred_dif = ((oto_size)^2 - new_oto_size6$mean))

fishdat6a <- fishdat6 %>% filter(pred_dif>0)
fishdat6b <- fishdat6 %>% filter(pred_dif<0)

high6 <- dplyr::summarise(fishdat6a, avg= mean(pred_dif))
low6 <- dplyr::summarise(fishdat6b, avg= mean(pred_dif))
print(high6)
print(low6)

q6 <- as.vector(quantile(fishdat6$pred_dif))


fishdat6 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()

#Graph of observed against predicted values for each
fishdat6 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat6 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik6 <- extract_log_lik(six_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo6 <- loo(lik6, save_psis = TRUE)
waic6 <- waic(lik6)
print(loo6)
plot(loo6)

loo::compare(loo4,loo6)


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

#For the original dataset
system.time(seven_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

#For the cut dataset
system.time(seven_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

#For the complete male dataset
system.time(seven_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(seven_fit, pars=c("alpha", "eta","beta1", "beta2", "epsilon","lp__", "yhat[110]", "slope_after", "intercept_after"))


##The wrong breakpoint is being discovered from this model
##I will cut the data off after this earlier breakpoint and then re-run the model 
##to see if the model can find the chosen breakpoint

##Extract the summary from the fit to plot the predicted lines
summary7 <- summary(seven_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  dplyr::select(variable, everything()) %>% 
  as_data_frame()

head(summary7)

new_oto_size7 <- summary7 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size7)

fishdat7 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size7$mean,
                            lower = new_oto_size7$`2.5%`,
                            upper = new_oto_size7$`97.5%`,
                            pred_dif = ((oto_size)^2 - new_oto_size7$mean))

fishdat7a <- fishdat7 %>% filter(pred_dif>0)
fishdat7b <- fishdat7 %>% filter(pred_dif<0)

high7 <- dplyr::summarise(fishdat7a, avg= mean(pred_dif))
low7 <- dplyr::summarise(fishdat7b, avg= mean(pred_dif))
print(high7)
print(low7)

q7 <- as.vector(quantile(fishdat7$pred_dif))

fishdat7 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour="steelblue", size=0.5) +
  geom_ribbon(aes(x=(prev)^2, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_classic()

#Graph of observed against predicted values for each
fishdat7 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat7 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()



##Want to extract the log-likelihood so that we can compare the models at some point
lik7 <- extract_log_lik(seven_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo7 <- loo(lik7, save_psis = TRUE)
waic7 <- waic(lik7)
print(loo7)
plot(loo7)

loo::compare(loo6,loo7)

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

#For the original dataset
system.time(eight_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the cut dataset
system.time(eight_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the complete male dataset
system.time(eight_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(eight_fit, pars=c("alpha[110]", "eta","beta1", "beta2", "epsilon","lp__", "sigma_alpha", "yhat[110]", "slope_after", "intercept_after[110]"))

##Extract the summary from the fit to plot the predicted lines
summary8 <- summary(eight_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  dplyr::select(variable, everything()) %>% 
  as_data_frame()

head(summary8)

new_oto_size8 <- summary8 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size8)

fishdat8 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size8$mean,
                                lower = new_oto_size8$`2.5%`,
                                upper = new_oto_size8$`97.5%`,
                                pred_dif = ((oto_size)^2 - new_oto_size8$mean))

fishdat8a <- fishdat8 %>% filter(pred_dif>0)
fishdat8b <- fishdat8 %>% filter(pred_dif<0)

high8 <- dplyr::summarise(fishdat8a, avg= mean(pred_dif))
low8 <- dplyr::summarise(fishdat8b, avg= mean(pred_dif))
print(high8)
print(low8)

q8 <- as.vector(quantile(fishdat8$pred_dif))

fishdat8 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour="steelblue", size=0.5) +
  geom_ribbon(aes(x=(prev)^2, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_classic()

#Graph of observed against predicted values for each
fishdat8 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat8 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik8 <- extract_log_lik(eight_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo8 <- loo(lik8, save_psis = TRUE)
waic8 <- waic(lik8)
print(loo8)
plot(loo8)

loo::compare(loo6,loo8)

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

#For the original dataset
system.time(nine_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the cut dataset
system.time(nine_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=11)))

#For the complete male dataset
system.time(nine_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(nine_fit, pars=c("alpha", "eta[110]", "mu_eta", "sigma_eta","beta1", "beta2", "epsilon","lp__", "yhat[110]","slope_after", "intercept_after[110]"))

##Extract the summary from the fit to plot the predicted lines
summary9 <- summary(nine_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  dplyr::select(variable, everything()) %>% 
  as_data_frame()

head(summary9)

new_oto_size9 <- summary9 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size9)

fishdat9 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size9$mean,
                                lower = new_oto_size9$`2.5%`,
                                upper = new_oto_size9$`97.5%`,
                                pred_dif = ((oto_size)^2 - new_oto_size9$mean))

fishdat9a <- fishdat9 %>% filter(pred_dif>0)
fishdat9b <- fishdat9 %>% filter(pred_dif<0)

high9 <- dplyr::summarise(fishdat9a, avg= mean(pred_dif))
low9 <- dplyr::summarise(fishdat9b, avg= mean(pred_dif))
print(high9)
print(low9)

q9 <- as.vector(quantile(fishdat9$pred_dif))

fishdat9 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour="steelblue", size=0.5) +
  geom_ribbon(aes(x=(prev)^2, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_classic()

#Graph of observed against predicted values for each
fishdat9 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat9 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik9 <- extract_log_lik(nine_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo9 <- loo(lik9, save_psis = TRUE)
waic9 <- waic(lik9)
print(loo9)
plot(loo9)

loo::compare(loo6,loo9)


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

#For the original dataset
system.time(ten_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the cut dataset
system.time(ten_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=11)))

#For the complete male dataset
system.time(ten_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(ten_fit, pars=c("alpha[110]", "eta[110]", "mu_eta", "sigma_eta","beta1", "beta2", "epsilon","lp__", "sigma_alpha", "yhat[110]", "slope_after", "intercept_after[110]"))

##Extract the summary from the fit to plot the predicted lines
summary10 <- summary(ten_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  dplyr::select(variable, everything()) %>% 
  as_data_frame()

head(summary10)

new_oto_size10 <- summary10 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size7)

fishdat10 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size10$mean,
                                lower = new_oto_size10$`2.5%`,
                                upper = new_oto_size10$`97.5%`,
                                pred_dif = ((oto_size)^2 - new_oto_size10$mean))

fishdat10a <- fishdat10 %>% filter(pred_dif>0)
fishdat10b <- fishdat10 %>% filter(pred_dif<0)

high10 <- dplyr::summarise(fishdat10a, avg= mean(pred_dif))
low10 <- dplyr::summarise(fishdat10b, avg= mean(pred_dif))
print(high10)
print(low10)

q10 <- as.vector(quantile(fishdat10$pred_dif))

fishdat10 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour="steelblue", size=0.5) +
  geom_ribbon(aes(x=(prev)^2, ymin = lower, ymax = upper), alpha = 0.25)+
  theme_classic()


#Graph of observed against predicted values for each
fishdat10 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat10 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik10 <- extract_log_lik(ten_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo10 <- loo(lik10, save_psis = TRUE)
waic10 <- waic(lik10)
print(loo10)
plot(loo10)

loo::compare(loo6,loo10)


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




###################################################################################
##Running the model. 11. Unhinged fixed threshold + intercept (varying slopes)####

rt <- stanc(file="The Stan models/Unhinged fixed threshold + intercept (varying slopes).stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

#For the original dataset
system.time(eleven_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=15, adapt_delta=0.8)))

#For the cut dataset
system.time(eleven_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

#For the complete male dataset
system.time(eleven_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))


##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(eleven_fit, pars=c("alpha1", "alpha2", "beta1", "beta2", "epsilon","lp__", "eta"))

##Extract the summary from the fit to plot the predicted lines
summary11 <- summary(eleven_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary11)

new_oto_size11 <- summary11 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size11)

fishdat11 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size11$mean,
                                lower = new_oto_size11$`2.5%`,
                                upper = new_oto_size11$`97.5%`,
                                pred_dif = ((oto_size)^2 - new_oto_size11$mean))

fishdat11a <- fishdat11 %>% filter(pred_dif>0)
fishdat11b <- fishdat11 %>% filter(pred_dif<0)

high11 <- dplyr::summarise(fishdat11a, avg= mean(pred_dif))
low11 <- dplyr::summarise(fishdat11b, avg= mean(pred_dif))
print(high11)
print(low11)

q11 <- as.vector(quantile(fishdat11$pred_dif))

fishdat11 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()


#Graph of observed against predicted values for each
fishdat11 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)

#Graph of oto-size against difference between predicted and observed data
fishdat11 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik11 <- extract_log_lik(eleven_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo11 <- loo(lik11, save_psis = TRUE)
waic11 <- waic(lik11)
print(loo11)
plot(loo11)

loo::compare(loo6,loo11)


#########################################################################################
##Running the model. 12. Unhinged fixed threshold, varying intercept (varying slopes)####

rt <- stanc(file="The Stan models/Unhinged fixed threshold, varying intercept (varying slopes).stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

#For the main dataset
system.time(twelve_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

#For the cut dataset
system.time(twelve_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=10, adapt_delta=0.8)))

#For the complete male dataset
system.time(twelve_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(twelve_fit, pars=c("alpha1[110]", "alpha2[110]", "eta","beta1", "beta2", "epsilon","lp__", "sigma_alpha", "yhat[110]"))

##Extract the summary from the fit to plot the predicted lines
summary12 <- summary(twelve_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary12)

new_oto_size12 <- summary12 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size12)

fishdat12 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size12$mean,
                                 lower = new_oto_size12$`2.5%`,
                                 upper = new_oto_size12$`97.5%`,
                                 pred_dif = ((oto_size)^2 - new_oto_size12$mean))

fishdat12a <- fishdat12 %>% filter(pred_dif>0)
fishdat12b <- fishdat12 %>% filter(pred_dif<0)

high12 <- dplyr::summarise(fishdat12a, avg= mean(pred_dif))
low12 <- dplyr::summarise(fishdat12b, avg= mean(pred_dif))
print(high12)
print(low12)

q12 <- as.vector(quantile(fishdat12$pred_dif))

fishdat12 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()

#Graph of observed against predicted values for each
fishdat12 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)

#Graph of oto-size against difference between predicted and observed data
fishdat12 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik12 <- extract_log_lik(twelve_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo12 <- loo(lik12, save_psis = TRUE)
waic12 <- waic(lik12)
print(loo12)
plot(loo12)

loo::compare(loo6,loo12)

#########################################################################################
##Running the model. 13. Unhinged fixed intercept varying threshold (varying slopes)####

rt <- stanc(file="The Stan models/Unhinged fixed intercept varying threshold (varying slopes).stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

#For the original dataset
system.time(thirteen_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the cut dataset
system.time(thirteen_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=11)))

#For the complete male dataset
system.time(thirteen_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(thirteen_fit, pars=c("alpha1", "alpha2", "eta[110]", "mu_eta", "sigma_eta","beta1", "beta2", "epsilon","lp__", "yhat[110]"))

##Extract the summary from the fit to plot the predicted lines
summary13 <- summary(thirteen_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary13)

new_oto_size13 <- summary13 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size13)

fishdat13 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size13$mean,
                                 lower = new_oto_size13$`2.5%`,
                                 upper = new_oto_size13$`97.5%`,
                                 pred_dif = ((oto_size)^2 - new_oto_size13$mean))

fishdat13a <- fishdat13 %>% filter(pred_dif>0)
fishdat13b <- fishdat13 %>% filter(pred_dif<0)

high13 <- dplyr::summarise(fishdat13a, avg= mean(pred_dif))
low13 <- dplyr::summarise(fishdat13b, avg= mean(pred_dif))
print(high13)
print(low13)

q13 <- as.vector(quantile(fishdat13$pred_dif))

fishdat13 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()


#Graph of observed against predicted values for each
fishdat13 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat13 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik13 <- extract_log_lik(thirteen_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo13 <- loo(lik13, save_psis = TRUE)
waic13 <- waic(lik13)
print(loo13)
plot(loo13)

loo::compare(loo6,loo13)

#######################################################################################
##Running the model. 14. Unhinged varying threshold and intercept (varying slopes)####

rt <- stanc(file="The Stan models/Unhinged varying threshold and intercept (varying slopes).stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

#For original dataset
system.time(fourteen_fit <- sampling(sm, data=fishdat, seed=1, iter=2000, chains=4, control = list(max_treedepth=10)))

#For the cut dataset
system.time(fourteen_fit <- sampling(sm, data=fishdat_cut, seed=1, iter=2000, chains=4, control = list(max_treedepth=11)))

#For the complete male dataset
system.time(fourteen_fit <- sampling(sm, data=fishdat_male, seed=1, iter=2000, chains=4, control = list(max_treedepth=11, adapt_delta=0.8)))

##Pairs plot to see any correlation amongst parameters that could lead to correlations that may lead to low E-BFMI 
##values (to understand which parameters should be altered)
pairs(fourteen_fit, pars=c("alpha1[110]", "alpha2[110]", "eta[110]", "mu_eta", "sigma_eta","beta1", "beta2", "epsilon","lp__", "sigma_alpha", "yhat[110]"))


##Extract the summary from the fit to plot the predicted lines
summary14 <- summary(fourteen_fit)$summary %>% 
  as.data.frame() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_data_frame()

head(summary14)

new_oto_size14 <- summary14 %>% 
  filter(str_detect(variable,'sim_oto_size') & !str_detect(variable,'log') & !str_detect(variable,'pp'))

head(new_oto_size14)

fishdat14 <- EBSm_cut %>% mutate(mean_oto_size = new_oto_size14$mean,
                                 lower = new_oto_size14$`2.5%`,
                                 upper = new_oto_size14$`97.5%`,
                                 pred_dif = ((oto_size)^2 - new_oto_size14$mean))

fishdat14a <- fishdat14 %>% filter(pred_dif>0)
fishdat14b <- fishdat14 %>% filter(pred_dif<0)

high14 <- dplyr::summarise(fishdat14a, avg= mean(pred_dif))
low14 <- dplyr::summarise(fishdat14b, avg= mean(pred_dif))
print(high14)
print(low14)

q14 <- as.vector(quantile(fishdat14$pred_dif))

fishdat14 %>% 
  ggplot() +
  geom_point(aes(x = (prev)^2, y = (oto_size)^2, colour=maturity), alpha=0.1) +
  geom_point(aes(x = (prev)^2, y = mean_oto_size), colour= "green", alpha=0.25, size=0.5 )+
  theme_classic()

#Graph of observed against predicted values for each
fishdat14 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = mean_oto_size), alpha=0.1)+
  theme_classic()

#Graph of oto-size against difference between predicted and observed data
fishdat14 %>% 
  ggplot() +
  geom_point(aes(x = (oto_size)^2, y = pred_dif), alpha=0.1)+
  geom_hline(yintercept=0)+
  theme_classic()

##Want to extract the log-likelihood so that we can compare the models at some point
lik14 <- extract_log_lik(fourteen_fit, parameter_name = "log_lik", merge_chains = TRUE)
loo14 <- loo(lik14, save_psis = TRUE)
waic14 <- waic(lik14)
print(loo14)
plot(loo14)

loo::compare(loo6,loo14)

loo_list <- list(loo1,loo2,loo3,loo4,loo5,loo6,loo7,loo8,loo9,loo10,loo11,loo12,loo13,loo14)
loo_model_weights(loo_list)
loo_model_weights(loo_list, method = "pseudobma")
loo_model_weights(loo_list, method = "pseudobma", BB=FALSE)



#Creating a box-plot to understand the difference between simulated data
#from each model and the 
quantiles <- cbind(fishdat1$pred_dif,fishdat2$pred_dif,fishdat3$pred_dif, fishdat4$pred_dif, fishdat5$pred_dif,fishdat6$pred_dif,
                   fishdat7$pred_dif,fishdat8$pred_dif,fishdat9$pred_dif,fishdat10$pred_dif,fishdat11$pred_dif,
                   fishdat12$pred_dif,fishdat13$pred_dif,fishdat14$pred_dif)
                   
              
quantiles <- as.data.frame(quantiles)

quantiles <- quantiles %>% tidyr::gather()

quantiles <- as.disc
quantiles %>% group_by(key)
quantiles %>% 
  ggplot() +
  geom_boxplot(aes(y=value, x=key))+
  geom_hline(yintercept =0, colour="red")+
  theme_classic()


#working with waics

waics <- c(
  waic1$estimates["elpd_waic", 1],
  waic2$estimates["elpd_waic", 1],
  waic3$estimates["elpd_waic", 1],
  waic4$estimates["elpd_waic", 1],
  waic5$estimates["elpd_waic", 1],
  waic6$estimates["elpd_waic", 1],
  waic7$estimates["elpd_waic", 1],
  waic8$estimates["elpd_waic", 1],
  waic9$estimates["elpd_waic", 1],
  waic10$estimates["elpd_waic", 1],
  waic11$estimates["elpd_waic", 1],
  waic12$estimates["elpd_waic", 1],
  waic13$estimates["elpd_waic", 1],
  waic14$estimates["elpd_waic", 1]
)

waics <- as.list(waics)
waic_wts <- waics /sum(waics)
round(waic_wts, 3) 
