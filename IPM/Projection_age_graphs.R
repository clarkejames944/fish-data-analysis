source("analysis/setup.R")

###################################################

## A function to add zone, temp and sex data to the summarised data
create_new_age_dataset <- function(IPM_summary, z, s){
  zone <- 1:nrow(IPM_summary)
  zone <- as.data.frame(zone)
  zone <- zone %>%  transmute(zone= c(z))
  sex <- 1:nrow(IPM_summary)
  sex <- as.data.frame(sex)
  sex <- sex %>%  transmute(sex= c(s))
  IPM_summary <- cbind(IPM_summary, zone)
  IPM_summary <- cbind(IPM_summary, sex)
}

################################################
###EBS summaries###########
IPM_name <- "EBSf_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
EBSf_res_summary <- load(file = f_loc)

EBSf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

EBSf_summary <- create_new_age_dataset(EBSf_summary, 'EBS', 'Female')


############################################
## male projection
IPM_name <- "EBSm_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
EBSm_res_summary <- load(file = f_loc)

EBSm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

EBSm_summary <- create_new_age_dataset(EBSm_summary, 'EBS', 'Male')


################################################
###ETAS summaries###########
IPM_name <- "ETASf_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
ETASf_res_summary <- load(file = f_loc)

ETASf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

ETASf_summary <- create_new_age_dataset(ETASf_summary, 'ETAS', 'Female')


##################################################################
## male projection
IPM_name <- "ETASm_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
ETASm_res_summary <- load(file = f_loc)

ETASm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

ETASm_summary <- create_new_age_dataset(ETASm_summary, 'ETAS', 'Male')


################################################
###WTAS summaries###########
IPM_name <- "WTASf_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
WTASf_res_summary <- load(file = f_loc)

WTASf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

WTASf_summary <- create_new_age_dataset(WTASf_summary, 'WTAS', 'Female')

##################################################################
## male projection
IPM_name <- "WTASm_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
WTASm_res_summary <- load(file = f_loc)

WTASm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

WTASm_summary <- create_new_age_dataset(WTASm_summary, 'WTAS', 'Male')


################################################
###NSW summaries###########
IPM_name <- "NSWf_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
NSWf_res_summary <- load(file = f_loc)

NSWf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

NSWf_summary <- create_new_age_dataset(NSWf_summary, 'NSW', 'Female')


##################################################################
## male projection
IPM_name <- "NSWm_age_at_maturity_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
NSWm_res_summary <- load(file = f_loc)

NSWm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_p   = mean(mu_p),
    q10_mu_p = quantile(mu_p, probs = 0.10), # for the 80% interval
    q75_mu_p = quantile(mu_p, probs = 0.75), # for the 50% interval
    q25_mu_p = quantile(mu_p, probs = 0.25), # for the 50% interval
    q90_mu_p = quantile(mu_p, probs = 0.90), # for the 80% interval
  )

NSWm_summary <- create_new_age_dataset(NSWm_summary, 'NSW', 'Male')


################################################
## Join the summary data
#####
## Join the datasets together
onejoin <- full_join(EBSf_summary, ETASf_summary)
twojoin <- full_join(onejoin, WTASf_summary)
full_f_summary <- full_join(twojoin, NSWf_summary)

## Now for the male
onejoin_m <- full_join(EBSm_summary, ETASm_summary)
twojoin_m <- full_join(onejoin_m, WTASm_summary)
full_m_summary <- full_join(twojoin_m, NSWm_summary)

## Join the two together
full_curr_summary <- full_join(full_f_summary, full_m_summary)

#################################################

ggplot(full_curr_summary, aes(x=a, y=mean_mu_p, linetype=sex))+
  geom_line()+
  geom_ribbon(aes(ymin=q25_mu_p, ymax=q75_mu_p, fill=sex), alpha=0.5)+
  facet_grid(~zone)+
  ylab("Proportion of individuals post-threshold") + xlab("Age (years)")+
  labs(fill= "Sex", linetype="Sex")+
  scale_fill_manual(values = c("royalblue3", "darkorange1"))+
  theme(axis.title = element_text(size=13),
         axis.text= element_text(size=11))
