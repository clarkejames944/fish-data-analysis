source("analysis/setup.R")

###################################################

## A function to add zone, temp and sex data to the summarised data
create_new_dataset <- function(IPM_summary, z, s, t){
  zone <- 1:nrow(IPM_summary)
  zone <- as.data.frame(zone)
  zone <- zone %>%  transmute(zone= c(z))
  sex <- 1:nrow(IPM_summary)
  sex <- as.data.frame(sex)
  sex <- sex %>%  transmute(sex= c(s))
  temp <- 1:nrow(IPM_summary)
  temp <- as.data.frame(temp)
  temp <- temp %>%  transmute(temp= c(t))
  IPM_summary <- cbind(IPM_summary, zone)
  IPM_summary <- cbind(IPM_summary, sex)
  IPM_summary <- cbind(IPM_summary, temp)
}

################################################
###EBS summaries###########
IPM_name <- "EBSf_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
EBSf_res_summary <- load(file = f_loc)

EBSf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

EBSf_summary <- create_new_dataset(EBSf_summary, 'EBS', 'Female', 'Current')

#########################################################
## For increased temperatures
IPM_name <- "EBSf+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
EBSf_temp_res_summary <- load(file = f_loc)

EBSf_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

EBSf_temp_summary <- create_new_dataset(EBSf_temp_summary, 'EBS', 'Female', 'Increased')


############################################
## male projection
IPM_name <- "EBSm_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
EBSm_res_summary <- load(file = f_loc)

EBSm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

EBSm_summary <- create_new_dataset(EBSm_summary, 'EBS', 'Male', 'Current')

##################################################################
## For increased temperatures
IPM_name <- "EBSm+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
EBSm_temp_res_summary <- load(file = f_loc)

EBSm_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

EBSm_temp_summary <- create_new_dataset(EBSm_temp_summary, 'EBS', 'Male', 'Increased')

################################################
###ETAS summaries###########
IPM_name <- "ETASf_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
ETASf_res_summary <- load(file = f_loc)

ETASf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

ETASf_summary <- create_new_dataset(ETASf_summary, 'ETAS', 'Female', 'Current')

##################################################################
## For increased temperatures
IPM_name <- "ETASf+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
ETASf_temp_res_summary <- load(file = f_loc)

ETASf_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

ETASf_temp_summary <- create_new_dataset(ETASf_temp_summary, 'ETAS', 'Female', 'Increased')

##################################################################
## male projection
IPM_name <- "ETASm_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
ETASm_res_summary <- load(file = f_loc)

ETASm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

ETASm_summary <- create_new_dataset(ETASm_summary, 'ETAS', 'Male', 'Current')

##################################################################
## For increased temperatures
IPM_name <- "ETASm+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
ETASm_temp_res_summary <- load(file = f_loc)

ETASm_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

ETASm_temp_summary <- create_new_dataset(ETASm_temp_summary, 'ETAS', 'Male', 'Increased')

################################################
###WTAS summaries###########
IPM_name <- "WTASf_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
WTASf_res_summary <- load(file = f_loc)

WTASf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

WTASf_summary <- create_new_dataset(WTASf_summary, 'WTAS', 'Female', 'Current')

##################################################################
## For increased temperatures
IPM_name <- "WTASf+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
WTASf_temp_res_summary <- load(file = f_loc)

WTASf_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

WTASf_temp_summary <- create_new_dataset(WTASf_temp_summary, 'WTAS', 'Female', 'Increased')

##################################################################
## male projection
IPM_name <- "WTASm_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
WTASm_res_summary <- load(file = f_loc)

WTASm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

WTASm_summary <- create_new_dataset(WTASm_summary, 'WTAS', 'Male', 'Current')

##################################################################
## For increased temperatures
IPM_name <- "WTASm+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
WTASm_temp_res_summary <- load(file = f_loc)

WTASm_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

WTASm_temp_summary <- create_new_dataset(WTASm_temp_summary, 'WTAS', 'Male', 'Increased')

################################################
###NSW summaries###########
IPM_name <- "NSWf_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
NSWf_res_summary <- load(file = f_loc)

NSWf_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

NSWf_summary <- create_new_dataset(NSWf_summary, 'NSW', 'Female', 'Current')

##################################################################

## For increased temperatures
IPM_name <- "NSWf+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
NSWf_temp_res_summary <- load(file = f_loc)

NSWf_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

NSWf_temp_summary <- create_new_dataset(NSWf_temp_summary, 'NSW', 'Female', 'Increased')

##################################################################
## male projection
IPM_name <- "NSWm_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
NSWm_res_summary <- load(file = f_loc)

NSWm_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

NSWm_summary <- create_new_dataset(NSWm_summary, 'NSW', 'Male', 'Current')

##################################################################
## For increased temperatures
IPM_name <- "NSWm+1_iteration_data"

f_loc <- sub("XX", replacement = IPM_name, x = "IPM_summaries/XX.rda")
NSWm_temp_res_summary <- load(file = f_loc)

NSWm_temp_summary <- res_summary %>% 
  group_by(a) %>% 
  summarise(
    mean_mu_z   = mean(mu_z),
    q10_mu_z = quantile(mu_z, probs = 0.10), # for the 80% interval
    q75_mu_z = quantile(mu_z, probs = 0.75), # for the 50% interval
    q25_mu_z = quantile(mu_z, probs = 0.25), # for the 50% interval
    q90_mu_z = quantile(mu_z, probs = 0.90), # for the 80% interval
  )

NSWm_temp_summary <- create_new_dataset(NSWm_temp_summary, 'NSW', 'Male', 'Increased')

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

##temp_join female
temp_onejoin <- full_join(EBSf_temp_summary, ETASf_temp_summary)
temp_twojoin <- full_join(temp_onejoin, WTASf_temp_summary)
temp_full_f_summary <- full_join(temp_twojoin, NSWf_temp_summary)

#temp_join male
temp_onejoin_m <- full_join(EBSm_temp_summary, ETASm_temp_summary)
temp_twojoin_m <- full_join(temp_onejoin_m, WTASm_temp_summary)
temp_full_m_summary <- full_join(temp_twojoin_m, NSWm_temp_summary)

###join the temp sexes together
full_temp_summary <- full_join(temp_full_f_summary, temp_full_m_summary)

### join it all
full_summary <- full_join(full_curr_summary, full_temp_summary)




######################################################################
## Make the plots

raw_size_at_age <- fishdat_cut %>% group_by(Age, zone, sex) %>% 
                      summarise(mean=mean((z0)^2),
                                se= sd((z0)^2)/sqrt(length(z0))
                                )
raw_size_at_age <- raw_size_at_age %>% transmute(a=Age, mean_mu_z=mean, se=se, sex=sex)

raw_size_at_age <- mutate(raw_size_at_age, sex= car::recode(sex, "'F' = 'Female'"))
raw_size_at_age <- mutate(raw_size_at_age, sex= car::recode(sex, "'M' = 'Male'"))

ggplot(full_curr_summary, aes(x=a, y=mean_mu_z, linetype=zone))+
  geom_line()+
  geom_ribbon(aes(ymin=q10_mu_z, ymax=q90_mu_z, fill=zone), alpha=0.35)+
  facet_wrap(~sex)+
  ylab("Mean otolith size (cm)") + xlab("Age (years)")+
  geom_point(data = raw_size_at_age, aes(color=zone), position = position_dodge(width = 0.5))+
  geom_errorbar(data = raw_size_at_age, aes(ymin=mean_mu_z-se, ymax=mean_mu_z+se, color=zone), position = position_dodge(width = 0.5))+
  theme_classic()



## plot of size at age separated by temperature and sex
full_summary %>% 
  ggplot(aes(x=a, y=mean_mu_z, linetype=sex))+
  geom_line(aes(color=temp))+
  geom_ribbon(aes(ymin=q10_mu_z, ymax=q90_mu_z, fill=temp), alpha=0.35)+
  facet_wrap(~zone)+
  ylab("Mean otolith size (cm)") + xlab("Age (years)")+
  geom_point(data = raw_size_at_age, position = position_dodge(width = 0.5))+
  geom_errorbar(data = raw_size_at_age, aes(ymin=mean_mu_z-se, ymax=mean_mu_z+se), position = position_dodge(width = 0.5))+
  theme_classic()


