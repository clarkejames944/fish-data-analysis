###############
####Time for some Stan####

#Setting up
#Libraries required
library(dplyr)
library(ggplot2)
library(rstan)
library(lattice)

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
set.seed(123)

##Specifying the data required for this model
fishmod <- lm(oto_size ~ prev, EBSm)
x <-model.matrix(fishmod)
fishdat <- list(N_EBSm=nrow(EBSm),
                Ngroups = length(unique(EBSm$FishID)),
                oto_size=EBSm$oto_size,
                Age= as.numeric(EBSm$Age),
                fishID=as.numeric(EBSm$FishID),
                prev=EBSm$prev
)

rt <- stanc(file="Stan code.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

View(EBSm)


