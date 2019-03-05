#Fish data analysis
#11/12/2018

#Setting up-----------
#Libraries required
library(dplyr)
library(ggplot2)
library(mgcv)
library(voxel)
library(MASS)
library(car)

#Data required
fish <- read.csv("https://raw.githubusercontent.com/clarkejames944/fish-data-analysis/master/otoliths%20(working)/data_derived/data_otolith_complete.csv")

#See the Data sleuthing notes Rmarkdown file for some general comments about 
#the datasets and their constituent variables

#explore the data a bit
glimpse(fish)
summary(fish)

str(fish)
summarise_all(fish, funs(mean))

#make a few plots of variables of interest- start with the effect of temperature on some variables
ggplot(fish, aes(x=bottomtemp1, y=Increment))+
  geom_point(size=1)

ggplot(fish, aes(x=bottomtemp1, y=log_incr))+
  geom_point(size=1)

ggplot(fish, aes(x=bottomtemp1, y=growth))+
  geom_point(size=1)

#Shows that as the temperature increases so does the 
#population-wide average growth


#make a few plots of variables of interest- have a look at age in years now

ggplot(fish, aes(x=Age, y=growth))+
  geom_point(size=1)

ggplot(fish, aes(x=AdjAge, y=radius))+
  geom_point(size=1)

ggplot(fish, aes(x=Age, y=Increment))+
  geom_point(size=1)

ggplot(fish, aes(x=Age, y=log_incr))+
  geom_point(size=1)

ggplot(fish, aes(x=AdjAge, y=floorlength))+
  geom_point(size=1)

#Look at time of year of capture- must group by a particular age before doing this one

ggplot(fish, aes(x=month, y=growth))+
  geom_point(size=1)

ggplot(fish, aes(x=month, y=radius))+
  geom_point(size=1)

ggplot(fish, aes(x=month, y=Increment))+
  geom_point(size=1)

ggplot(fish, aes(x=month, y=log_incr))+
  geom_point(size=1)

#Have a look at gender

ggplot(fish, aes(x=sex, y=growth))+
  geom_boxplot()

ggplot(fish, aes(x=sex, y=radius))+
  geom_boxplot()

ggplot(fish, aes(x=sex, y=Increment))+
  geom_boxplot()

ggplot(fish, aes(x=sex, y=log_incr))+
  geom_boxplot()

ggplot(fish, aes(x=sex, y=AdjAge))+
  geom_boxplot()

#Effect of area caught

ggplot(fish, aes(x=area, y=growth))+
  geom_boxplot()

ggplot(fish, aes(x=area, y=radius))+
  geom_boxplot()

ggplot(fish, aes(x=area, y=Increment))+
  geom_boxplot()

ggplot(fish, aes(x=area, y=log_incr))+
  geom_boxplot()

ggplot(fish, aes(x=area, y=floorlength))+
  geom_boxplot()

#Effect of zone caught

ggplot(fish, aes(x=zone, y=growth))+
  geom_boxplot()

ggplot(fish, aes(x=zone, y=radius))+
  geom_boxplot()

ggplot(fish, aes(x=zone, y=Increment))+
  geom_boxplot()

ggplot(fish, aes(x=zone, y=log_incr))+
  geom_boxplot()

ggplot(fish, aes(x=zone, y=floorlength))+
  geom_boxplot()

#plotting fish length against otolith size variables

ggplot(fish, aes(x=floorlength, y=growth))+
  geom_point(size=1)
#The above one doesn't match up in terms one is continuous across 
#a fishes life and the other is discrete

ggplot(fish, aes(x=floorlength, y=radius))+
  geom_point(size=1)

ggplot(fish, aes(x=floorlength, y=Increment))+
  geom_point(size=1)
#Again doesn't match up

ggplot(fish, aes(x=floorlength, y=log_incr))+
  geom_point(size=1)
#Again doesn't match

#Fishing gear effects

ggplot(fish, aes(x=gear, y=growth))+
  geom_boxplot()

ggplot(fish, aes(x=gear, y=radius))+
  geom_boxplot()

ggplot(fish, aes(x=gear, y=Increment))+
  geom_boxplot()

ggplot(fish, aes(x=gear, y=log_incr))+
  geom_boxplot()


#######################################
###Fitting some GAMs to the data####

#Setting up-----------
#Libraries required
library(dplyr)
library(ggplot2)
library(mgcv)
library(voxel)
library(MASS)
library(car)
library(ggfortify)
library(scales)
library(devtools)

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


##Now all is set up to do a GAM
##Oto_size against previous year
?gam


#create subsets


EBSm <- not_1 %>% filter(sex=='M', zone=='EBS')
EBSf <- not_1 %>% filter(sex=='F', zone=='EBS')
ETASm <- not_1 %>% filter(sex=='M', zone=='ETAS')
ETASf <- not_1 %>% filter(sex=='F', zone=='ETAS')
NSWm <- not_1 %>% filter(sex=='M', zone=='NSW')
NSWf <- not_1 %>% filter(sex=='F', zone=='NSW')
WTASm <- not_1 %>% filter(sex=='M', zone=='WTAS')
WTASf <- not_1 %>% filter(sex=='F', zone=='WTAS')

##I will start with the EBS male data and create a GAM from this
EBSm_size_gam <- gam(EBSm$oto_size~s(EBSm$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(EBSm_size_gam)

summary(EBSm_size_gam)

#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function

#plot the graph for EBSm
EBSm <- EBSm %>% mutate(preds=predict(EBSm_size_gam))
ggplot(EBSm,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("EBSm")+
  theme_classic()



##GAM between size and extracted residuals
##extract residuals
EBSm_res_orig<- residuals(EBSm_size_gam)
#absolute values only
EBSm_res_size <- (EBSm_res_orig)^2
EBSm <- EBSm %>% mutate(res_oto_size=EBSm_res_size)

#The GAM of size against the extracted residuals
EBSm_res_gam <- gam(EBSm_res_size~s(EBSm$prev))
gam.check(EBSm_res_gam)
summary(EBSm_res_gam)

EBSm <- EBSm %>% mutate(res_preds=predict(EBSm_res_gam))
ggplot(EBSm,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("EBSm")+
  theme_classic()

##Same again for EBS females
EBSf_size_gam <- gam(EBSf$oto_size~s(EBSf$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(EBSf_size_gam)

summary(EBSf_size_gam)

#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function

EBSf <- EBSf %>% mutate(preds=predict(EBSf_size_gam))
ggplot(EBSf,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("EBSf")+
  theme_classic()



##GAM between size and extracted residuals
##extract residuals
EBSf_res_orig<- residuals(EBSf_size_gam)
#absolute values only
EBSf_res_size <- (EBSf_res_orig)^2
EBSf <- EBSf %>% mutate(res_oto_size=EBSf_res_size)

#The GAM of size against the extracted residuals
EBSf_res_gam <- gam(EBSf_res_size~s(EBSf$prev))
gam.check(EBSf_res_gam)
summary(EBSf_res_gam)

EBSf <- EBSf %>% mutate(res_preds=predict(EBSf_res_gam))
ggplot(EBSf,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("EBSf")+
  theme_classic()

##now for ETAS males
ETASm_size_gam <- gam(ETASm$oto_size~s(ETASm$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(ETASm_size_gam)

summary(ETASm_size_gam)

#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function
ETASm <- ETASm %>% mutate(preds=predict(ETASm_size_gam))

ggplot(ETASm,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("ETASm")+
  theme_classic()


##GAM between size and extracted residuals
##extract residuals
ETASm_res_orig<- residuals(ETASm_size_gam)
#absolute values only
ETASm_res_size <- (ETASm_res_orig)^2
ETASm <- ETASm %>% mutate(res_oto_size=ETASm_res_size)

#The GAM of size against the extracted residuals
ETASm_res_gam <- gam(ETASm_res_size~s(ETASm$prev))
gam.check(ETASm_res_gam)
summary(ETASm_res_gam)

ETASm <- ETASm %>% mutate(res_preds=predict(ETASm_res_gam))
ggplot(ETASm,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("ETASm")+
  theme_classic()


##now for ETAS females
ETASf_size_gam <- gam(ETASf$oto_size~s(ETASf$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(ETASf_size_gam)

summary(ETASf_size_gam)

#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function
ETASf <- ETASf %>% mutate(preds=predict(ETASf_size_gam))

ggplot(ETASf,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("ETASf")+
  theme_classic()



##GAM between size and extracted residuals
##extract residuals
ETASf_res_orig<- residuals(ETASf_size_gam)
#absolute values only
ETASf_res_size <- (ETASf_res_orig)^2
ETASf <- ETASf %>% mutate(res_oto_size=ETASf_res_size)

#The GAM of size against the extracted residuals
ETASf_res_gam <- gam(ETASf_res_size~s(ETASf$prev))
gam.check(ETASf_res_gam)
summary(ETASf_res_gam)

ETASf <- ETASf %>% mutate(res_preds=predict(ETASf_res_gam))
ggplot(ETASf,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("ETASf")+
  theme_classic()

##now for NSW males
NSWm_size_gam <- gam(NSWm$oto_size~s(NSWm$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(NSWm_size_gam)

summary(NSWm_size_gam)

#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function
NSWm <- NSWm %>% mutate(preds=predict(NSWm_size_gam))
ggplot(NSWm,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("NSWm")+
  theme_classic()


##GAM between size and extracted residuals
##extract residuals
NSWm_res_orig<- residuals(NSWm_size_gam)
#absolute values only
NSWm_res_size <-(NSWm_res_orig)^2
NSWm <- NSWm %>% mutate(res_oto_size=NSWm_res_size)

#The GAM of size against the extracted residuals
NSWm_res_gam <- gam(NSWm_res_size~s(NSWm$prev))
gam.check(NSWm_res_gam)
summary(NSWm_res_gam)

NSWm <- NSWm %>% mutate(res_preds=predict(NSWm_res_gam))
ggplot(NSWm,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("NSWm")+
  theme_classic()


##now for NSW females
NSWf_size_gam <- gam(NSWf$oto_size~s(NSWf$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(NSWf_size_gam)

summary(NSWf_size_gam)

NSWf <- NSWf %>% mutate(preds=predict(NSWf_size_gam))
ggplot(NSWf,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("NSWf")+
  theme_classic()


##GAM between size and extracted residuals
##extract residuals
NSWf_res_orig<- residuals(NSWf_size_gam)
#absolute values only
NSWf_res_size <- (NSWf_res_orig)^2
NSWf <- NSWf %>% mutate(res_oto_size=NSWf_res_size)

#The GAM of size against the extracted residuals
NSWf_res_gam <- gam(NSWf_res_size~s(NSWf$prev))
gam.check(NSWf_res_gam)
summary(NSWf_res_gam)

NSWf <- NSWf %>% mutate(res_preds=predict(NSWf_res_gam))
ggplot(NSWf,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("NSWf")+
  theme_classic()

##now for WTAS males
WTASm_size_gam <- gam(WTASm$oto_size~s(WTASm$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(WTASm_size_gam)

summary(WTASm_size_gam)

#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function
WTASm <- WTASm %>% mutate(preds=predict(WTASm_size_gam))
ggplot(WTASm,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("WTASm")+
  theme_classic()


##GAM between size and extracted residuals
##extract residuals
WTASm_res_orig<- residuals(WTASm_size_gam)
#absolute values only
WTASm_res_size <- (WTASm_res_orig)^2
WTASm <- WTASm %>% mutate(res_oto_size=WTASm_res_size)

#The GAM of size against the extracted residuals
WTASm_res_gam <- gam(WTASm_res_size~s(WTASm$prev))
gam.check(WTASm_res_gam)
summary(WTASm_res_gam)

WTASm <- WTASm %>% mutate(res_preds=predict(WTASm_res_gam))
ggplot(WTASm,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("WTASm")+
  theme_classic()


##now for WTAS females
WTASf_size_gam <- gam(WTASf$oto_size~s(WTASf$prev))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(WTASf_size_gam)

summary(WTASf_size_gam)

#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function

WTASf <- WTASf %>% mutate(preds=predict(WTASf_size_gam))
ggplot(WTASf,aes(x=prev, y=oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("WTASf")+
  theme_classic()


##GAM between size and extracted residuals
##extract residuals
WTASf_res_orig<- residuals(WTASf_size_gam)
#absolute values only
WTASf_res_size <- (WTASf_res_orig)^2
WTASf <- WTASf %>% mutate(res_oto_size=WTASf_res_size)

#The GAM of size against the extracted residuals
WTASf_res_gam <- gam(WTASf_res_size~s(WTASf$prev))
gam.check(WTASf_res_gam)
summary(WTASf_res_gam)

WTASf <- WTASf %>% mutate(res_preds=predict(WTASf_res_gam))
ggplot(WTASf,aes(x=prev, y=res_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("WTASf")+
  theme_classic()








##########Box-Cox Transformations#######
##Now we can move on to the transformations of these
#This shows that as size increases the residuals are decreasing 
#We would like to remove this effect as much as possible
#So we will try box-cox transformations

##Box-cox transformations of the otolith size data


##Try a manual method for Box-Cox transformation
##Vary lambda across a range of values

oto_box1=(NSWm$oto_size^-2 -1)/-2
oto_box2=(NSWm$oto_size^-1.5 -1)/-1.5
oto_box3=(NSWm$oto_size^-1 -1)/-1
oto_box4=(NSWm$oto_size^-0.5 -1)/-0.5
oto_box6=(NSWm$oto_size^0.5 -1)/0.5
oto_box7=(NSWm$oto_size^1 -1)/1
oto_box8=(NSWm$oto_size^1.5 -1)/1.5
oto_box9=(NSWm$oto_size^2 -1)/2
#Make a GAM of this new boxcox data
#Fit this over a range of values see how it changes


#create the box-cox GAM

NSWm_box_size_gam1 <- gam(oto_box1~s(NSWm$prev))
gam.check(NSWm_box_size_gam1)
summary(NSWm_box_size_gam1)
NSWm <- NSWm %>% mutate(box_preds1=predict(NSWm_box_size_gam1))
ggplot(NSWm,aes(x=prev, y=oto_box1))+
  geom_point(size=1)+
  geom_line(aes(x=oto_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm1")+
  theme_classic()

NSWm_box_size_gam2 <- gam(oto_box2~s(NSWm$prev))
gam.check(NSWm_box_size_gam2)
summary(NSWm_box_size_gam2)
NSWm <- NSWm %>% mutate(box_preds2=predict(NSWm_box_size_gam2))
ggplot(NSWm,aes(x=prev, y=oto_box2))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds2), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm2")+
  theme_classic()

NSWm_box_size_gam3 <- gam(oto_box3~s(NSWm$prev))
gam.check(NSWm_box_size_gam3)
summary(NSWm_box_size_gam3)
NSWm <- NSWm %>% mutate(box_preds3=predict(NSWm_box_size_gam3))
ggplot(NSWm,aes(x=prev, y=oto_box3))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds3), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm3")+
  theme_classic()

NSWm_box_size_gam4 <- gam(oto_box4~s(NSWm$prev))
gam.check(NSWm_box_size_gam4)
summary(NSWm_box_size_gam4)
NSWm <- NSWm %>% mutate(box_preds4=predict(NSWm_box_size_gam4))
ggplot(NSWm,aes(x=prev, y=oto_box4))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds4), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm4")+
  theme_classic()

NSWm_box_size_gam5 <- gam(oto_box5~s(NSWm$prev))
gam.check(NSWm_box_size_gam5)
summary(NSWm_box_size_gam5)
NSWm <- NSWm %>% mutate(box_preds5=predict(NSWm_box_size_gam5))
ggplot(NSWm,aes(x=prev, y=oto_box5))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds5), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm5")+
  theme_classic()

NSWm_box_size_gam6 <- gam(oto_box6~s(NSWm$prev))
gam.check(NSWm_box_size_gam6)
summary(NSWm_box_size_gam6)
NSWm <- NSWm %>% mutate(box_preds6=predict(NSWm_box_size_gam6))
ggplot(NSWm,aes(x=prev, y=oto_box6))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds6), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm6")+
  theme_classic()

NSWm_box_size_gam7 <- gam(oto_box7~s(NSWm$prev))
gam.check(NSWm_box_size_gam7)
summary(NSWm_box_size_gam7)
NSWm <- NSWm %>% mutate(box_preds7=predict(NSWm_box_size_gam7))
ggplot(NSWm,aes(x=prev, y=oto_box7))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds7), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm7")+
  theme_classic()

NSWm_box_size_gam8 <- gam(oto_box8~s(NSWm$prev))
gam.check(NSWm_box_size_gam8)
summary(NSWm_box_size_gam8)
NSWm <- NSWm %>% mutate(box_preds8=predict(NSWm_box_size_gam8))
ggplot(NSWm,aes(x=prev, y=oto_box8))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds8), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm8")+
  theme_classic()

NSWm_box_size_gam9 <- gam(oto_box9~s(NSWm$prev))
gam.check(NSWm_box_size_gam9)
summary(NSWm_box_size_gam9)
NSWm <- NSWm %>% mutate(box_preds9=predict(NSWm_box_size_gam9))
ggplot(NSWm,aes(x=prev, y=oto_box9))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds9), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm9")+
  theme_classic()

#Extract the residuals from these

NSWm_res_box_size1 <- residuals(NSWm_box_size_gam1)
NSWm_res_box_size1 <- (NSWm_res_box_size1)^2
NSWm <- NSWm %>% mutate(res_box_size1=NSWm_res_box_size1)

NSWm_res_box_size2 <- residuals(NSWm_box_size_gam2)
NSWm_res_box_size2 <- (NSWm_res_box_size2)^2
NSWm <- NSWm %>% mutate(res_box_size2=NSWm_res_box_size2)

NSWm_res_box_size3 <- residuals(NSWm_box_size_gam3)
NSWm_res_box_size3 <- (NSWm_res_box_size3)^2
NSWm <- NSWm %>% mutate(res_box_size3=NSWm_res_box_size3)

NSWm_res_box_size4 <- residuals(NSWm_box_size_gam4)
NSWm_res_box_size4 <- (NSWm_res_box_size4)^2
NSWm <- NSWm %>% mutate(res_box_size4=NSWm_res_box_size4)

NSWm_res_box_size6 <- residuals(NSWm_box_size_gam6)
NSWm_res_box_size6 <- (NSWm_res_box_size6)^2
NSWm <- NSWm %>% mutate(res_box_size6=NSWm_res_box_size6)

NSWm_res_box_size7 <- residuals(NSWm_box_size_gam7)
NSWm_res_box_size7 <- (NSWm_res_box_size7)^2
NSWm <- NSWm %>% mutate(res_box_size7=NSWm_res_box_size7)

NSWm_res_box_size8 <- residuals(NSWm_box_size_gam8)
NSWm_res_box_size8 <- (NSWm_res_box_size8)^2
NSWm <- NSWm %>% mutate(res_box_size8=NSWm_res_box_size8)

NSWm_res_box_size9 <- residuals(NSWm_box_size_gam9)
NSWm_res_box_size9 <- (NSWm_res_box_size9)^2
NSWm <- NSWm %>% mutate(res_box_size9=NSWm_res_box_size9)

#GAM of the residuals

NSWm_res_box_gam1 <- gam(NSWm$res_box_size1~s(NSWm$prev))
gam.check(NSWm_res_box_gam1)
summary(NSWm_res_box_gam1)
NSWm <- NSWm %>% mutate(res_box_preds1=predict(NSWm_res_box_gam1))
ggplot(NSWm,aes(x=prev, y=res_box_size1))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds1), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm1")+
  theme_classic()

NSWm_res_box_gam2 <- gam(NSWm$res_box_size2~s(NSWm$prev))
gam.check(NSWm_res_box_gam2)
summary(NSWm_res_box_gam2)
NSWm <- NSWm %>% mutate(res_box_preds2=predict(NSWm_res_box_gam2))
ggplot(NSWm,aes(x=prev, y=res_box_size2))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds2), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm2")+
  theme_classic()

NSWm_res_box_gam3 <- gam(NSWm$res_box_size3~s(NSWm$prev))
gam.check(NSWm_res_box_gam3)
summary(NSWm_res_box_gam3)
NSWm <- NSWm %>% mutate(res_box_preds3=predict(NSWm_res_box_gam3))
ggplot(NSWm,aes(x=prev, y=res_box_size3))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds3), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm3")+
  theme_classic()

NSWm_res_box_gam4 <- gam(NSWm$res_box_size4~s(NSWm$prev))
gam.check(NSWm_res_box_gam4)
summary(NSWm_res_box_gam4)
NSWm <- NSWm %>% mutate(res_box_preds4=predict(NSWm_res_box_gam4))
ggplot(NSWm,aes(x=prev, y=res_box_size4))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds4), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm4")+
  theme_classic()

NSWm_res_box_gam6 <- gam(NSWm$res_box_size6~s(NSWm$prev))
gam.check(NSWm_res_box_gam6)
summary(NSWm_res_box_gam6)
NSWm <- NSWm %>% mutate(res_box_preds6=predict(NSWm_res_box_gam6))
ggplot(NSWm,aes(x=prev, y=res_box_size6))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds6), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm1")+
  theme_classic()

NSWm_res_box_gam7 <- gam(NSWm$res_box_size7~s(NSWm$prev))
gam.check(NSWm_res_box_gam7)
summary(NSWm_res_box_gam7)
NSWm <- NSWm %>% mutate(res_box_preds7=predict(NSWm_res_box_gam7))
ggplot(NSWm,aes(x=prev, y=res_box_size7))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds7), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm7")+
  theme_classic()

NSWm_res_box_gam8 <- gam(NSWm$res_box_size8~s(NSWm$prev))
gam.check(NSWm_res_box_gam8)
summary(NSWm_res_box_gam8)
NSWm <- NSWm %>% mutate(res_box_preds8=predict(NSWm_res_box_gam8))
ggplot(NSWm,aes(x=prev, y=res_box_size8))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds8), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm8")+
  theme_classic()

NSWm_res_box_gam9 <- gam(NSWm$res_box_size9~s(NSWm$prev))
gam.check(NSWm_res_box_gam9)
summary(NSWm_res_box_gam9)
NSWm <- NSWm %>% mutate(res_box_preds9=predict(NSWm_res_box_gam9))
ggplot(NSWm,aes(x=prev, y=res_box_size9))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds9), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWm8")+
  theme_classic()


##On a graph altogether these lines for each value of lambda are:

ggplot(NSWm, aes(x=prev))+
  geom_line(aes(y=res_box_preds1), colour="steelblue")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds1), label="lambda=-2")+
  geom_line(aes(y=res_box_preds2), colour="green")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds2), label="lambda=-1.5")+
  geom_line(aes(y=res_box_preds3), colour="red")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds3), label="lambda=-1")+
  geom_line(aes(y=res_box_preds4), colour="orange")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds4), label="lambda=-0.5")+
  geom_line(aes(y=res_box_preds6), colour="aquamarine")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds6), label="lambda=0.5")+
  geom_line(aes(y=res_box_preds7), colour="mediumorchid1")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds7), label="lambda=1")+
  geom_line(aes(y=res_box_preds8), colour="purple")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds8), label="lambda=1.5")+
  geom_line(aes(y=res_box_preds9), colour="cyan")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds9), label="lambda=2")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox NSWm")+
  theme_classic()

##Seems like lambda=1-2 is best
##Keep these ones in the plot
ggplot(NSWm, aes(x=prev))+
  geom_line(aes(y=res_box_preds6), colour="aquamarine")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds6), label="lambda=0.5")+
  geom_line(aes(y=res_box_preds7), colour="mediumorchid1")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds7), label="lambda=1")+
  geom_line(aes(y=res_box_preds8), colour="purple")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds8), label="lambda=1.5")+
  geom_line(aes(y=res_box_preds9), colour="cyan")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds9), label="lambda=2")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox NSWm")+
  theme_classic()

##Somewhere in between 1.5 and 2 would be good for lambda
oto_boxA=(NSWm$oto_size^1.75 -1)/1.75

NSWm_box_size_gamA <- gam(oto_boxA~s(NSWm$prev))
gam.check(NSWm_box_size_gamA)
summary(NSWm_box_size_gamA)
NSWm <- NSWm %>% mutate(box_predsA=predict(NSWm_box_size_gamA))
ggplot(NSWm,aes(x=prev, y=oto_boxA))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_predsA), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmA")+
  theme_classic()

NSWm_res_box_sizeA <- residuals(NSWm_box_size_gamA)
NSWm_res_box_sizeA <- (NSWm_res_box_sizeA)^2
NSWm <- NSWm %>% mutate(res_box_sizeA=NSWm_res_box_sizeA)

NSWm_res_box_gamA <- gam(NSWm$res_box_sizeA~s(NSWm$prev))
gam.check(NSWm_res_box_gamA)
summary(NSWm_res_box_gamA)
NSWm <- NSWm %>% mutate(res_box_predsA=predict(NSWm_res_box_gamA))
ggplot(NSWm,aes(x=prev, y=res_box_sizeA))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_predsA), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmA")+
  theme_classic()

##compare this with the others
ggplot(NSWm, aes(x=prev))+
  geom_line(aes(y=res_box_preds8), colour="purple")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds8), label="lambda=1.5")+
  geom_line(aes(y=res_box_preds9), colour="cyan")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds9), label="lambda=2")+
  geom_line(aes(y=res_box_predsA), colour="olivedrab")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_predsA), label="lambda=1.75")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox NSWm")+
  theme_classic()

#Between 2 and 1.75
#Try 1.85
oto_boxB=(NSWm$oto_size^1.85 -1)/1.85

NSWm_box_size_gamB <- gam(oto_boxB~s(NSWm$prev))
gam.check(NSWm_box_size_gamB)
summary(NSWm_box_size_gamB)
NSWm <- NSWm %>% mutate(box_predsB=predict(NSWm_box_size_gamB))
ggplot(NSWm,aes(x=prev, y=oto_boxB))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_predsB), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmB")+
  coord_fixed(ratio = 1)+
  theme_classic()

NSWm_res_box_sizeB <- residuals(NSWm_box_size_gamB)
NSWm_res_box_sizeB <- (NSWm_res_box_sizeB)^2
NSWm <- NSWm %>% mutate(res_box_sizeB=NSWm_res_box_sizeB)

NSWm_res_box_gamB <- gam(NSWm$res_box_sizeB~s(NSWm$prev, k=4))
gam.check(NSWm_res_box_gamB)
summary(NSWm_res_box_gamB)
NSWm <- NSWm %>% mutate(res_box_predsB=predict(NSWm_res_box_gamB))
ggplot(NSWm,aes(x=prev, y=res_box_sizeB))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_predsB), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmB")+
  coord_fixed(ratio = 30)+
  theme_classic()

ggplot(NSWm, aes(x=prev))+
  geom_line(aes(y=res_box_preds8), colour="purple")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds8), label="lambda=1.5")+
  geom_line(aes(y=res_box_preds9), colour="cyan")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds9), label="lambda=2")+
  geom_line(aes(y=res_box_predsA), colour="olivedrab")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_predsA), label="lambda=1.75")+
  geom_line(aes(y=res_box_predsB), colour="gray")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_predsB), label="lambda=1.85")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox NSWm")+
  theme_classic()

##1.85 is the best value of lambda so far
##optimum will be somewhere between 1.85 and 2 for lambda

#try 1.9
oto_boxC=(NSWm$oto_size^1.9 -1)/1.9

NSWm_box_size_gamC <- gam(oto_boxC~s(NSWm$prev))
gam.check(NSWm_box_size_gamC)
summary(NSWm_box_size_gamC)
NSWm <- NSWm %>% mutate(box_predsC=predict(NSWm_box_size_gamC))
ggplot(NSWm,aes(x=prev, y=oto_boxC))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_predsC), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmC")+
  theme_classic()

NSWm_res_box_sizeC <- residuals(NSWm_box_size_gamC)
NSWm_res_box_sizeC <- (NSWm_res_box_sizeC)^2
NSWm <- NSWm %>% mutate(res_box_sizeC=NSWm_res_box_sizeC)

NSWm_res_box_gamC <- gam(NSWm$res_box_sizeC~s(NSWm$prev))
gam.check(NSWm_res_box_gamC)
summary(NSWm_res_box_gamC)
NSWm <- NSWm %>% mutate(res_box_predsC=predict(NSWm_res_box_gamC))
ggplot(NSWm,aes(x=prev, y=res_box_sizeC))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_predsC), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmC")+
  theme_classic()

ggplot(NSWm, aes(x=prev))+
  geom_line(aes(y=res_box_preds8), colour="purple")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds8), label="lambda=1.5")+
  geom_line(aes(y=res_box_preds9), colour="cyan")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_preds9), label="lambda=2")+
  geom_line(aes(y=res_box_predsA), colour="olivedrab")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_predsA), label="lambda=1.75")+
  geom_line(aes(y=res_box_predsB), colour="gray")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_predsB), label="lambda=1.85")+
  geom_line(aes(y=res_box_predsC), colour="steelblue")+
  geom_text(data=subset(NSWm, prev=='0.5'), aes(y=res_box_predsC), label="lambda=1.9")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox NSWm")+
  theme_classic()
  
  ##Answer is somewhere between 1.85 and 1.9
  #I'll stick with lambda=1.85

#So for NSWm with lamnda=1.85
oto_boxB=(NSWm$oto_size^1.85 -1)/1.85
NSWm_box_prev = (NSWm$prev^1.85 -1)/1.85

NSWm_box_size_gamB <- gam(oto_boxB~s(NSWm$prev))
gam.check(NSWm_box_size_gamB)
summary(NSWm_box_size_gamB)
NSWm <- NSWm %>% mutate(box_predsB=predict(NSWm_box_size_gamB))
ggplot(NSWm,aes(x=prev, y=oto_boxB))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_predsB), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmB")+
  coord_fixed(ratio = 1)+
  theme_classic()

NSWm_res_box_orig <- residuals(NSWm_box_size_gamB)
NSWm_res_box_sizeB <- (NSWm_res_box_orig)^2
NSWm <- NSWm %>% mutate(res_box_sizeB=NSWm_res_box_sizeB)

NSWm_res_box_gamB <- gam(NSWm$res_box_sizeB~s(NSWm$prev))
gam.check(NSWm_res_box_gamB)
summary(NSWm_res_box_gamB)
NSWm <- NSWm %>% mutate(res_box_predsB=predict(NSWm_res_box_gamB))
ggplot(NSWm,aes(x=prev, y=res_box_sizeB))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_predsB), size=1.3, colour="steelblue")+
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmB")+
  coord_fixed(ratio = 30)+
  theme_classic()

##transforming x variable
NSWm_box_size_gam1B <- gam(oto_boxB~s(NSWm_box_prev))
gam.check(NSWm_box_size_gam1B)
summary(NSWm_box_size_gam1B)
NSWm <- NSWm %>% mutate(box_preds1B=predict(NSWm_box_size_gam1B))
ggplot(NSWm,aes(x=NSWm_box_prev, y=oto_boxB))+
  geom_point(size=1)+
  geom_line(aes(x=NSWm_box_prev, y=box_preds1B), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmB")+
  coord_fixed(ratio = 1)+
  theme_classic()

NSWm_res_box_orig1 <- residuals(NSWm_box_size_gam1B)
NSWm_res_box_size1B <- (NSWm_res_box_orig1)^2
NSWm <- NSWm %>% mutate(res_box_size1B=NSWm_res_box_size1B)

NSWm_res_box_gam1B <- gam(NSWm$res_box_size1B~s(NSWm_box_prev))
gam.check(NSWm_res_box_gam1B)
summary(NSWm_res_box_gam1B)
NSWm <- NSWm %>% mutate(res_box_preds1B=predict(NSWm_res_box_gam1B))
ggplot(NSWm,aes(x=NSWm_box_prev, y=res_box_size1B))+
  geom_point(size=1)+
  geom_line(aes(x=NSWm_box_prev, y=res_box_preds1B), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("Box Cox s")+
  ggtitle("Manual Box-Cox NSWmB")+
  coord_fixed(ratio = 30)+
  theme_classic()

  
  ##For NSWf
  NSWf_oto_box = (NSWf$oto_size^1.85 -1)/1.85
  NSWf_box_prev = (NSWf$prev^1.85 -1)/1.85
  
  NSWf_box_size_gam <- gam(NSWf_oto_box~s(NSWf$prev))
  gam.check(NSWf_box_size_gam)
  summary(NSWf_box_size_gam)
  NSWf <- NSWf %>% mutate(box_preds=predict(NSWf_box_size_gam))
  ggplot(NSWf,aes(x=prev, y=NSWf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox NSWf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  NSWf_res_box_orig <- residuals(NSWf_box_size_gam)
  NSWf_res_box_size <- (NSWf_res_box_orig)^2
  NSWf <- NSWf %>% mutate(res_box_size=NSWf_res_box_size)
  
  NSWf_res_box_gam <- gam(NSWf$res_box_size~s(NSWf_box_prev))
  gam.check(NSWf_res_box_gam)
  summary(NSWf_res_box_gam)
  NSWf <- NSWf %>% mutate(res_box_preds=predict(NSWf_res_box_gam))
  ggplot(NSWf,aes(x=prev, y=res_box_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox NSWf")+
    coord_fixed(ratio = 30)+
    theme_classic()
 
  #try transforming the x variable as well
  NSWf_box_size_gam1 <- gam(NSWf_oto_box~s(NSWf_box_prev))
  gam.check(NSWf_box_size_gam1)
  summary(NSWf_box_size_gam1)
  NSWf <- NSWf %>% mutate(box_preds1=predict(NSWf_box_size_gam1))
  ggplot(NSWf,aes(x=NSWf_box_prev, y=NSWf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=NSWf_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox NSWf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  NSWf_res_box_orig1 <- residuals(NSWf_box_size_gam1)
  NSWf_res_box_size1 <- (NSWf_res_box_orig1)^2
  NSWf <- NSWf %>% mutate(res_box_size1=NSWf_res_box_size1)
  
  NSWf_res_box_gam1 <- gam(NSWf$res_box_size1~s(NSWf_box_prev))
  gam.check(NSWf_res_box_gam1)
  summary(NSWf_res_box_gam1)
  NSWf <- NSWf %>% mutate(res_box_preds1=predict(NSWf_res_box_gam1))
  ggplot(NSWf,aes(x=NSWf_box_prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=NSWf_box_prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox NSWf")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##For EBSm
  EBSm_oto_box = (EBSm$oto_size^1.85 -1)/1.85
  EBSm_box_prev = (EBSm$prev^1.85 -1)/1.85
  
  EBSm_box_size_gam <- gam(EBSm_oto_box~s(EBSm$prev))
  gam.check(EBSm_box_size_gam)
  summary(EBSm_box_size_gam)
  EBSm <- EBSm %>% mutate(box_preds=predict(EBSm_box_size_gam))
  ggplot(EBSm,aes(x=prev, y=EBSm_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSm")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  EBSm_res_box_orig <- residuals(EBSm_box_size_gam)
  EBSm_res_box_size <- (EBSm_res_box_orig)^2
  EBSm <- EBSm %>% mutate(res_box_size=EBSm_res_box_size)
  
  EBSm_res_box_gam <- gam(EBSm$res_box_size~s(EBSm$prev))
  gam.check(EBSm_res_box_gam)
  summary(EBSm_res_box_gam)
  EBSm <- EBSm %>% mutate(res_box_preds=predict(EBSm_res_box_gam))
  ggplot(EBSm,aes(x=prev, y=res_box_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSm")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##Transforming x variable
  EBSm_oto_box = (EBSm$oto_size^1.85 -1)/1.85
  EBSm_box_prev = (EBSm$prev^1.85 -1)/1.85
  
  EBSm_box_size_gam1 <- gam(EBSm_oto_box~s(EBSm_box_prev))
  gam.check(EBSm_box_size_gam1)
  summary(EBSm_box_size_gam1)
  EBSm <- EBSm %>% mutate(box_preds1=predict(EBSm_box_size_gam1))
  ggplot(EBSm,aes(x=EBSm_box_prev, y=EBSm_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=EBSm_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSm")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  EBSm_res_box_orig1 <- residuals(EBSm_box_size_gam1)
  EBSm_res_box_size1 <- (EBSm_res_box_orig1)^2
  EBSm <- EBSm %>% mutate(res_box_size1=EBSm_res_box_size1)
  
  EBSm_res_box_gam1 <- gam(EBSm$res_box_size1~s(EBSm_box_prev))
  gam.check(EBSm_res_box_gam1)
  summary(EBSm_res_box_gam1)
  EBSm <- EBSm %>% mutate(res_box_preds1=predict(EBSm_res_box_gam1))
  ggplot(EBSm,aes(x=EBSm_box_prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=EBSm_box_prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSm")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##For EBSf
  
  EBSf_oto_box = (EBSf$oto_size^1.85 -1)/1.85
  EBSf_box_prev = (EBSf$prev^1.85 -1)/1.85
  
  EBSf_box_size_gam <- gam(EBSf_oto_box~s(EBSf$prev))
  gam.check(EBSf_box_size_gam)
  summary(EBSf_box_size_gam)
  EBSf <- EBSf %>% mutate(box_preds=predict(EBSf_box_size_gam))
  ggplot(EBSf,aes(x=prev, y=EBSf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  EBSf_res_box_orig <- residuals(EBSf_box_size_gam)
  EBSf_res_box_size <- (EBSf_res_box_orig)^2
  EBSf <- EBSf %>% mutate(res_box_size=EBSf_res_box_size)
  
  EBSf_res_box_gam <- gam(EBSf$res_box_size~s(EBSf$prev))
  gam.check(EBSf_res_box_gam)
  summary(EBSf_res_box_gam)
  EBSf <- EBSf %>% mutate(res_box_preds=predict(EBSf_res_box_gam))
  ggplot(EBSf,aes(x=prev, y=res_box_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSf")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##Transforming x variable
  EBSf_box_size_gam1 <- gam(EBSf_oto_box~s(EBSf_box_prev))
  gam.check(EBSf_box_size_gam1)
  summary(EBSf_box_size_gam1)
  EBSf <- EBSf %>% mutate(box_preds1=predict(EBSf_box_size_gam1))
  ggplot(EBSf,aes(x=EBSf_box_prev, y=EBSf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=EBSf_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  EBSf_res_box_orig1 <- residuals(EBSf_box_size_gam1)
  EBSf_res_box_size1 <- (EBSf_res_box_orig1)^2
  EBSf <- EBSf %>% mutate(res_box_size1=EBSf_res_box_size1)
  
  EBSf_res_box_gam1 <- gam(EBSf$res_box_size1~s(EBSf_box_prev))
  gam.check(EBSf_res_box_gam1)
  summary(EBSf_res_box_gam1)
  EBSf <- EBSf %>% mutate(res_box_preds1=predict(EBSf_res_box_gam1))
  ggplot(EBSf,aes(x=EBSf_box_prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=EBSf_box_prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox EBSf")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  
  ##For ETASm
  
  ETASm_oto_box = (ETASm$oto_size^1.85 -1)/1.85
  ETASm_box_prev = (ETASm$prev^1.85 -1)/1.85
  
  ETASm_box_size_gam <- gam(ETASm_oto_box~s(ETASm$prev))
  gam.check(ETASm_box_size_gam)
  summary(ETASm_box_size_gam)
  ETASm <- ETASm %>% mutate(box_preds=predict(ETASm_box_size_gam))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  ETASm_res_box_orig <- residuals(ETASm_box_size_gam)
  ETASm_res_box_size <- (ETASm_res_box_orig)^2
  ETASm <- ETASm %>% mutate(res_box_size=ETASm_res_box_size)
  
  ETASm_res_box_gam <- gam(ETASm$res_box_size~s(ETASm$prev))
  gam.check(ETASm_res_box_gam)
  summary(ETASm_res_box_gam)
  ETASm <- ETASm %>% mutate(res_box_preds=predict(ETASm_res_box_gam))
  ggplot(ETASm,aes(x=prev, y=res_box_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##Transforming x variable
  ETASm_box_size_gam1 <- gam(ETASm_oto_box~s(ETASm_box_prev))
  gam.check(ETASm_box_size_gam1)
  summary(ETASm_box_size_gam1)
  ETASm <- ETASm %>% mutate(box_preds1=predict(ETASm_box_size_gam1))
  ggplot(ETASm,aes(x=ETASm_box_prev, y=ETASm_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=ETASm_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  ETASm_res_box_orig1 <- residuals(ETASm_box_size_gam1)
  ETASm_res_box_size1 <- (ETASm_res_box_orig1)^2
  ETASm <- ETASm %>% mutate(res_box_size1=ETASm_res_box_size1)
  
  ETASm_res_box_gam1 <- gam(ETASm$res_box_size1~s(ETASm_box_prev))
  gam.check(ETASm_res_box_gam1)
  summary(ETASm_res_box_gam1)
  ETASm <- ETASm %>% mutate(res_box_preds1=predict(ETASm_res_box_gam1))
  ggplot(ETASm,aes(x=ETASm_box_prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=ETASm_box_prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##For ETASf
  
  ETASf_oto_box = (ETASf$oto_size^1.85 -1)/1.85
  ETASf_box_prev = (ETASf$prev^1.85 -1)/1.85
  
  ETASf_box_size_gam <- gam(ETASf_oto_box~s(ETASf$prev))
  gam.check(ETASf_box_size_gam)
  summary(ETASf_box_size_gam)
  ETASf <- ETASf %>% mutate(box_preds=predict(ETASf_box_size_gam))
  ggplot(ETASf,aes(x=prev, y=ETASf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  ETASf_res_box_orig <- residuals(ETASf_box_size_gam)
  ETASf_res_box_size <- (ETASf_res_box_orig)^2
  ETASf <- ETASf %>% mutate(res_box_size=ETASf_res_box_size)
  
  ETASf_res_box_gam <- gam(ETASf$res_box_size~s(ETASf$prev))
  gam.check(ETASf_res_box_gam)
  summary(ETASf_res_box_gam)
  ETASf <- ETASf %>% mutate(res_box_preds=predict(ETASf_res_box_gam))
  ggplot(ETASf,aes(x=prev, y=res_box_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASf")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##Transforming x variable
  
  ETASf_box_size_gam1 <- gam(ETASf_oto_box~s(ETASf_box_prev))
  gam.check(ETASf_box_size_gam1)
  summary(ETASf_box_size_gam1)
  ETASf <- ETASf %>% mutate(box_preds1=predict(ETASf_box_size_gam1))
  ggplot(ETASf,aes(x=ETASf_box_prev, y=ETASf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=ETASf_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  ETASf_res_box_orig1 <- residuals(ETASf_box_size_gam1)
  ETASf_res_box_size1 <- (ETASf_res_box_orig1)^2
  ETASf <- ETASf %>% mutate(res_box_size1=ETASf_res_box_size1)
  
  ETASf_res_box_gam1 <- gam(ETASf$res_box_size1~s(ETASf_box_prev))
  gam.check(ETASf_res_box_gam1)
  summary(ETASf_res_box_gam1)
  ETASf <- ETASf %>% mutate(res_box_preds1=predict(ETASf_res_box_gam1))
  ggplot(ETASf,aes(x=ETASf_box_prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=ETASf_box_prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASf")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##For WTASm
  
  WTASm_oto_box = (WTASm$oto_size^1.85 -1)/1.85
  WTASm_box_prev = (WTASm$prev^1.85 -1)/1.85
  
  WTASm_box_size_gam <- gam(WTASm_oto_box~s(WTASm$prev))
  gam.check(WTASm_box_size_gam)
  summary(WTASm_box_size_gam)
  WTASm <- WTASm %>% mutate(box_preds=predict(WTASm_box_size_gam))
  ggplot(WTASm,aes(x=prev, y=WTASm_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASm")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  WTASm_res_box_orig <- residuals(WTASm_box_size_gam)
  WTASm_res_box_size <- (WTASm_res_box_orig)^2
  WTASm <- WTASm %>% mutate(res_box_size=WTASm_res_box_size)
  
  WTASm_res_box_gam <- gam(WTASm$res_box_size~s(WTASm$prev))
  gam.check(WTASm_res_box_gam)
  summary(WTASm_res_box_gam)
  WTASm <- WTASm %>% mutate(res_box_preds=predict(WTASm_res_box_gam))
  ggplot(WTASm,aes(x=prev, y=res_box_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASm")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##Transforming x variable
  WTASm_box_size_gam1 <- gam(WTASm_oto_box~s(WTASm_box_prev))
  gam.check(WTASm_box_size_gam1)
  summary(WTASm_box_size_gam1)
  WTASm <- WTASm %>% mutate(box_preds1=predict(WTASm_box_size_gam1))
  ggplot(WTASm,aes(x=WTASm_box_prev, y=WTASm_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=WTASm_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASm")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  WTASm_res_box_orig1 <- residuals(WTASm_box_size_gam1)
  WTASm_res_box_size1 <- (WTASm_res_box_orig1)^2
  WTASm <- WTASm %>% mutate(res_box_size1=WTASm_res_box_size1)
  
  WTASm_res_box_gam1 <- gam(WTASm$res_box_size1~s(WTASm_box_prev))
  gam.check(WTASm_res_box_gam1)
  summary(WTASm_res_box_gam1)
  WTASm <- WTASm %>% mutate(res_box_preds1=predict(WTASm_res_box_gam1))
  ggplot(WTASm,aes(x=WTASm_box_prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=WTASm_box_prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASm")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##For WTASf
  WTASf_oto_box = (WTASf$oto_size^1.85 -1)/1.85
  WTASf_box_prev = (WTASf$prev^1.85 -1)/1.85
  
  WTASf_box_size_gam <- gam(WTASf_oto_box~s(WTASf$prev))
  gam.check(WTASf_box_size_gam)
  summary(WTASf_box_size_gam)
  WTASf <- WTASf %>% mutate(box_preds=predict(WTASf_box_size_gam))
  ggplot(WTASf,aes(x=prev, y=WTASf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  WTASf_res_box_orig <- residuals(WTASf_box_size_gam)
  WTASf_res_box_size <- (WTASf_res_box_orig)^2
  WTASf <- WTASf %>% mutate(res_box_size=WTASf_res_box_size)
  
  WTASf_res_box_gam <- gam(WTASf$res_box_size~s(WTASf$prev))
  gam.check(WTASf_res_box_gam)
  summary(WTASf_res_box_gam)
  WTASf <- WTASf %>% mutate(res_box_preds=predict(WTASf_res_box_gam))
  ggplot(WTASf,aes(x=prev, y=res_box_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASf")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##Transforming x variable
  
  WTASf_box_size_gam1 <- gam(WTASf_oto_box~s(WTASf_box_prev))
  gam.check(WTASf_box_size_gam1)
  summary(WTASf_box_size_gam1)
  WTASf <- WTASf %>% mutate(box_preds1=predict(WTASf_box_size_gam1))
  ggplot(WTASf,aes(x=WTASf_box_prev, y=WTASf_oto_box))+
    geom_point(size=1)+
    geom_line(aes(x=WTASf_box_prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASf")+
    coord_fixed(ratio = 1)+
    theme_classic()
  
  WTASf_res_box_orig1 <- residuals(WTASf_box_size_gam1)
  WTASf_res_box_size1 <- (WTASf_res_box_orig1)^2
  WTASf <- WTASf %>% mutate(res_box_size1=WTASf_res_box_size1)
  
  WTASf_res_box_gam1 <- gam(WTASf$res_box_size1~s(WTASf_box_prev))
  gam.check(WTASf_res_box_gam1)
  summary(WTASf_res_box_gam1)
  WTASf <- WTASf %>% mutate(res_box_preds1=predict(WTASf_res_box_gam1))
  ggplot(WTASf,aes(x=WTASf_box_prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=WTASf_box_prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox WTASf")+
    coord_fixed(ratio = 30)+
    theme_classic()
  
  ##I'll have a go at improving the ETASm transformation
  
  ETASm_oto_box1=(ETASm$oto_size^-2 -1)/-2
  ETASm_oto_box2=(ETASm$oto_size^-1.5 -1)/-1.5
  ETASm_oto_box3=(ETASm$oto_size^-1 -1)/-1
  ETASm_oto_box4=(ETASm$oto_size^-0.5 -1)/-0.5
  ETASm_oto_box6=(ETASm$oto_size^0.5 -1)/0.5
  ETASm_oto_box7=(ETASm$oto_size^1 -1)/1
  ETASm_oto_box8=(ETASm$oto_size^1.5 -1)/1.5
  ETASm_oto_box9=(ETASm$oto_size^2 -1)/2
  #Make a GAM of this new boxcox data
  #Fit this over a range of values see how it changes
  
  #create the box-cox GAM for each value of lambda
 
  ETASm_box_size_gam1 <- gam(ETASm_oto_box1~s(ETASm$prev))
  gam.check(ETASm_box_size_gam1)
  summary(ETASm_box_size_gam1)
  ETASm <- ETASm %>% mutate(box_preds1=predict(ETASm_box_size_gam1))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box1))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm1")+
    theme_classic()
  
  ETASm_box_size_gam2 <- gam(ETASm_oto_box2~s(ETASm$prev))
  gam.check(ETASm_box_size_gam2)
  summary(ETASm_box_size_gam2)
  ETASm <- ETASm %>% mutate(box_preds2=predict(ETASm_box_size_gam2))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box2))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds2), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm2")+
    theme_classic()
  
  ETASm_box_size_gam3 <- gam(ETASm_oto_box3~s(ETASm$prev))
  gam.check(ETASm_box_size_gam3)
  summary(ETASm_box_size_gam3)
  ETASm <- ETASm %>% mutate(box_preds3=predict(ETASm_box_size_gam3))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box3))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds3), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm3")+
    theme_classic()
  
  ETASm_box_size_gam4 <- gam(ETASm_oto_box4~s(ETASm$prev))
  gam.check(ETASm_box_size_gam4)
  summary(ETASm_box_size_gam4)
  ETASm <- ETASm %>% mutate(box_preds4=predict(ETASm_box_size_gam4))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box4))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds4), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm4")+
    theme_classic()
  
  ETASm_box_size_gam5 <- gam(ETASm_oto_box5~s(ETASm$prev))
  gam.check(ETASm_box_size_gam5)
  summary(ETASm_box_size_gam5)
  ETASm <- ETASm %>% mutate(box_preds5=predict(ETASm_box_size_gam5))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box5))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds5), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm5")+
    theme_classic()
  
  ETASm_box_size_gam6 <- gam(ETASm_oto_box6~s(ETASm$prev))
  gam.check(ETASm_box_size_gam6)
  summary(ETASm_box_size_gam6)
  ETASm <- ETASm %>% mutate(box_preds6=predict(ETASm_box_size_gam6))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box6))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds6), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm6")+
    theme_classic()
  
  ETASm_box_size_gam7 <- gam(ETASm_oto_box7~s(ETASm$prev))
 
  ETASm <- ETASm %>% mutate(box_preds7=predict(ETASm_box_size_gam7))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box7))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds7), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm7")+
    theme_classic()
  
  ETASm_box_size_gam8 <- gam(ETASm_oto_box8~s(ETASm$prev))
 
  ETASm <- ETASm %>% mutate(box_preds8=predict(ETASm_box_size_gam8))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box8))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds8), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm8")+
    theme_classic()
  
  ETASm_box_size_gam9 <- gam(ETASm_oto_box9~s(ETASm$prev))
  
  ETASm <- ETASm %>% mutate(box_preds9=predict(ETASm_box_size_gam9))
  ggplot(ETASm,aes(x=prev, y=ETASm_oto_box9))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=box_preds9), size=1.3, colour="steelblue")+
    ylab("Box Cox s'")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm9")+
    theme_classic()
  
  #Extract the residuals from these for each value of lambda

   ETASm_res_box_size1 <- residuals(ETASm_box_size_gam1)
  ETASm_res_box_size1 <- (ETASm_res_box_size1)^2
  ETASm <- ETASm %>% mutate(res_box_size1=ETASm_res_box_size1)
  
  ETASm_res_box_size2 <- residuals(ETASm_box_size_gam2)
  ETASm_res_box_size2 <- (ETASm_res_box_size2)^2
  ETASm <- ETASm %>% mutate(res_box_size2=ETASm_res_box_size2)
  
  ETASm_res_box_size3 <- residuals(ETASm_box_size_gam3)
  ETASm_res_box_size3 <- (ETASm_res_box_size3)^2
  ETASm <- ETASm %>% mutate(res_box_size3=ETASm_res_box_size3)
  
  ETASm_res_box_size4 <- residuals(ETASm_box_size_gam4)
  ETASm_res_box_size4 <- (ETASm_res_box_size4)^2
  ETASm <- ETASm %>% mutate(res_box_size4=ETASm_res_box_size4)
  
  ETASm_res_box_size6 <- residuals(ETASm_box_size_gam6)
  ETASm_res_box_size6 <- (ETASm_res_box_size6)^2
  ETASm <- ETASm %>% mutate(res_box_size6=ETASm_res_box_size6)
  
  ETASm_res_box_size7 <- residuals(ETASm_box_size_gam7)
  ETASm_res_box_size7 <- (ETASm_res_box_size7)^2
  ETASm <- ETASm %>% mutate(res_box_size7=ETASm_res_box_size7)
  
  ETASm_res_box_size8 <- residuals(ETASm_box_size_gam8)
  ETASm_res_box_size8 <- (ETASm_res_box_size8)^2
  ETASm <- ETASm %>% mutate(res_box_size8=ETASm_res_box_size8)
  
  ETASm_res_box_size9 <- residuals(ETASm_box_size_gam9)
  ETASm_res_box_size9 <- (ETASm_res_box_size9)^2
  ETASm <- ETASm %>% mutate(res_box_size9=ETASm_res_box_size9)
  
  #GAM of the residuals for each size of lambda
ETASm_res_box_gam1 <- gam(ETASm$res_box_size1~s(ETASm$prev))
  gam.check(ETASm_res_box_gam1)
  summary(ETASm_res_box_gam1)
  ETASm <- ETASm %>% mutate(res_box_preds1=predict(ETASm_res_box_gam1))
  ggplot(ETASm,aes(x=prev, y=res_box_size1))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds1), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm1")+
    theme_classic()
  
  ETASm_res_box_gam2 <- gam(ETASm$res_box_size2~s(ETASm$prev))
  gam.check(ETASm_res_box_gam2)
  summary(ETASm_res_box_gam2)
  ETASm <- ETASm %>% mutate(res_box_preds2=predict(ETASm_res_box_gam2))
  ggplot(ETASm,aes(x=prev, y=res_box_size2))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds2), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm2")+
    theme_classic()
  
  ETASm_res_box_gam3 <- gam(ETASm$res_box_size3~s(ETASm$prev))
  gam.check(ETASm_res_box_gam3)
  summary(ETASm_res_box_gam3)
  ETASm <- ETASm %>% mutate(res_box_preds3=predict(ETASm_res_box_gam3))
  ggplot(ETASm,aes(x=prev, y=res_box_size3))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds3), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm3")+
    theme_classic()
  
  ETASm_res_box_gam4 <- gam(ETASm$res_box_size4~s(ETASm$prev))
  gam.check(ETASm_res_box_gam4)
  summary(ETASm_res_box_gam4)
  ETASm <- ETASm %>% mutate(res_box_preds4=predict(ETASm_res_box_gam4))
  ggplot(ETASm,aes(x=prev, y=res_box_size4))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds4), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm4")+
    theme_classic()
  
  ETASm_res_box_gam6 <- gam(ETASm$res_box_size6~s(ETASm$prev))
  gam.check(ETASm_res_box_gam6)
  summary(ETASm_res_box_gam6)
  ETASm <- ETASm %>% mutate(res_box_preds6=predict(ETASm_res_box_gam6))
  ggplot(ETASm,aes(x=prev, y=res_box_size6))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds6), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm1")+
    theme_classic()
  
  ETASm_res_box_gam7 <- gam(ETASm$res_box_size7~s(ETASm$prev))
  gam.check(ETASm_res_box_gam7)
  summary(ETASm_res_box_gam7)
  ETASm <- ETASm %>% mutate(res_box_preds7=predict(ETASm_res_box_gam7))
  ggplot(ETASm,aes(x=prev, y=res_box_size7))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds7), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm7")+
    theme_classic()
  
  ETASm_res_box_gam8 <- gam(ETASm$res_box_size8~s(ETASm$prev))
  gam.check(ETASm_res_box_gam8)
  summary(ETASm_res_box_gam8)
  ETASm <- ETASm %>% mutate(res_box_preds8=predict(ETASm_res_box_gam8))
  ggplot(ETASm,aes(x=prev, y=res_box_size8))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds8), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm8")+
    theme_classic()
  
  ETASm_res_box_gam9 <- gam(ETASm$res_box_size9~s(ETASm$prev))
  gam.check(ETASm_res_box_gam9)
  summary(ETASm_res_box_gam9)
  ETASm <- ETASm %>% mutate(res_box_preds9=predict(ETASm_res_box_gam9))
  ggplot(ETASm,aes(x=prev, y=res_box_size9))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_box_preds9), size=1.3, colour="steelblue")+
    ylab("Box Cox r^2")+
    xlab("Box Cox s")+
    ggtitle("Manual Box-Cox ETASm9")+
    theme_classic()
  
  
  ##On a graph altogether these lines from each model are:

  ggplot(ETASm, aes(x=prev))+
    geom_line(aes(y=res_box_preds1), colour="steelblue")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds1), label="lambda=-2")+
    geom_line(aes(y=res_box_preds2), colour="green")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds2), label="lambda=-1.5")+
    geom_line(aes(y=res_box_preds3), colour="red")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds3), label="lambda=-1")+
    geom_line(aes(y=res_box_preds4), colour="orange")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds4), label="lambda=-0.5")+
    geom_line(aes(y=res_box_preds6), colour="aquamarine")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds6), label="lambda=0.5")+
    geom_line(aes(y=res_box_preds7), colour="mediumorchid1")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds7), label="lambda=1")+
    geom_line(aes(y=res_box_preds8), colour="purple")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds8), label="lambda=1.5")+
    geom_line(aes(y=res_box_preds9), colour="cyan")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds9), label="lambda=2")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("Manual Box-Cox ETASm")+
    theme_classic()
 
  ##Only the last three lines now
  ##checking between 1 and 2 for the lambda value
  ggplot(ETASm, aes(x=prev))+
    geom_line(aes(y=res_box_preds7), colour="mediumorchid1")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds7), label="lambda=1")+
    geom_line(aes(y=res_box_preds8), colour="orange")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds8), label="lambda=1.5")+
    geom_line(aes(y=res_box_preds9), colour="cyan")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds9), label="lambda=2")+
    geom_line(aes(y=res_box_preds), colour="steelblue")+
    geom_text(data=subset(ETASm, prev=='0.343'), aes(y=res_box_preds), label="lambda=1.85")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("Manual Box-Cox ETASm")+
    theme_classic()
  
  ##None of these seem to want to flatten 
  ##Still the one at 1.5 seems to be the best of the three
  ##Could probably stick with 1.85
  
  #####qq-plots####
  
  #carry out some qq plots for the residuals from each subset
  #NSWf and transformed NSWf
qqnorm(EBSm_res_orig)
  qqline(EBSm_res_orig)
  
  qqnorm(EBSm_res_box_orig)
  qqline(EBSm_res_box_orig)
  
  qqnorm(EBSm_res_box_orig1)
  qqline(EBSm_res_box_orig1)
  
    
  #EBSf and transformed EBSf
  qqnorm(EBSf_res_orig)
  qqline(EBSf_res_orig)
  
  qqnorm(EBSf_res_box_orig)
  qqline(EBSf_res_box_orig)
  
  qqnorm(EBSf_res_box_orig1)
  qqline(EBSf_res_box_orig1)
  

#ETASm and transformed ETASm
  qqnorm(ETASm_res_orig)
  qqline(ETASm_res_orig)
  
  qqnorm(ETASm_res_box_orig)
  qqline(ETASm_res_box_orig)
  
  qqnorm(ETASm_res_box_orig1)
  qqline(ETASm_res_box_orig1)
  
#ETASf and transformed ETASf

  qqnorm(ETASf_res_orig)
  qqline(ETASf_res_orig)
  
  qqnorm(ETASf_res_box_orig)
  qqline(ETASf_res_box_orig)
  
  qqnorm(ETASf_res_box_orig1)
  qqline(ETASf_res_box_orig1)
  

#NSWm and transformed NSWm
  qqnorm(NSWm_res_orig)
  qqline(NSWm_res_orig)
  
  qqnorm(NSWm_res_box_orig)
  qqline(NSWm_res_box_orig)
  
  qqnorm(NSWm_res_box_orig1)
  qqline(NSWm_res_box_orig1)
 

 #NSWf and transformed NSWf
  qqnorm(NSWf_res_orig)
  qqline(NSWf_res_orig)
  
  qqnorm(NSWf_res_box_orig)
  qqline(NSWf_res_box_orig)
  
  qqnorm(NSWf_res_box_orig1)
  qqline(NSWf_res_box_orig1)

  
#WTASm and transformed WTASm
  qqnorm(WTASm_res_orig)
  qqline(WTASm_res_orig)
  
  qqnorm(WTASm_res_box_orig)
  qqline(WTASm_res_box_orig)
  
  qqnorm(WTASm_res_box_orig1)
  qqline(WTASm_res_box_orig1)
  

#WTASf and transformed WTASf
  qqnorm(WTASf_res_orig)
  qqline(WTASf_res_orig)
  
  qqnorm(WTASf_res_box_orig)
  qqline(WTASf_res_box_orig)
  
  qqnorm(WTASf_res_box_orig1)
  qqline(WTASf_res_box_orig1)

##All of these q-qplots are right-skewed
#Therefore it may be time to specify another distribution for the model.
  
#############Time to try another model.#######
  
  #Setting up-----------
  #Libraries required
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(voxel)
  library(MASS)
  library(car)
  library(ggfortify)
  library(scales)
  library(devtools)
  
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
  
  
  ##Now all is set up to do a GAM
  ##Oto_size against previous year
  ?gam
  
  
  #create subsets
  
  
  EBSm <- not_1 %>% filter(sex=='M', zone=='EBS')
  EBSf <- not_1 %>% filter(sex=='F', zone=='EBS')
  ETASm <- not_1 %>% filter(sex=='M', zone=='ETAS')
  ETASf <- not_1 %>% filter(sex=='F', zone=='ETAS')
  NSWm <- not_1 %>% filter(sex=='M', zone=='NSW')
  NSWf <- not_1 %>% filter(sex=='F', zone=='NSW')
  WTASm <- not_1 %>% filter(sex=='M', zone=='WTAS')
  WTASf <- not_1 %>% filter(sex=='F', zone=='WTAS')

  ##Back to the raw data and I'll try specifying a gamma distribution first
  ##I'll see what happens with the EBSm data first
  
  
  
  ##I will start with the EBS male data and create a GAM from this
  EBSm_size_ggam <- gam(EBSm$oto_size~EBSm$prev+s(EBSm$prev, k=50, by=EBSm$maturity), family=Gamma(link="inverse"))
  gam.check(EBSm_size_ggam)
  summary(EBSm_size_ggam)
  #plot the graph for EBSm
  EBSm <- EBSm %>% mutate(gpreds=predict(EBSm_size_ggam, type="response"))
  ggplot(EBSm,aes(x=prev))+
    geom_point(aes(y=oto_size, colour=maturity), size=1, alpha=0.1)+
    geom_abline(intercept=0, slope= 1)+
    geom_point(aes(y=gpreds), size=0.25, colour="steelblue")+
    ylab("s'")+ xlab("s")+
   ggtitle("EBSm")+
    theme_classic()
  
  
  ##extract residuals
  EBSm_res_orig = residuals(EBSm_size_ggam, type="response")
  ###log the residuals
  EBSm_res_size = log(EBSm_res_orig^2)

  
  ggplot(EBSm,aes(x=gpreds, y=EBSm_res_size))+
    geom_point(aes(colour=maturity), size=1, alpha=0.1)+
    geom_smooth(method='gam',
                formula = y~s(x, k=100), se=FALSE)+
    geom_smooth(method='lm',
                se=FALSE, linetype=2)

 
   #GAM of residuals for normality    
  EBSm_res_ggam <- gam(EBSm$res_oto_size~s(EBSm$gpreds))
  gam.check(EBSm_res_ggam)
  summary(EBSm_res_ggam)
  
##Make some scaled residuals
  EBSm <- mutate(EBSm,
                       sc_resid = EBSm_res_orig / sqrt(exp(predict(EBSm_res_ggam))))
  
  ##Then make some qqplots (separting adult and juvenile)
  with(EBSm, car::qqp(sc_resid[maturity == "juvenile"])) 
  with(EBSm, car::qqp(sc_resid[maturity == "adult"])) 

  
  
  
  ##Same again for EBS females
  EBSf_size_gam <- gam(EBSf$oto_size~s(EBSf$prev, k=4))
  ##using gam.check helps to influence the choice of the value of k
  ##4 seems to be the best to use in this case
  gam.check(EBSf_size_gam)
  
  summary(EBSf_size_gam)
  
  #R-sq=0.978
  #Deviance explained = 97.8%
  #GCV=0.0018687
  #GCV is the measure of the degree of smoothness of the function
  
  EBSf <- EBSf %>% mutate(preds=predict(EBSf_size_gam))
  ggplot(EBSf,aes(x=prev, y=oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
    ylab("s'")+
    xlab("s")+
    ggtitle("EBSf")+
    theme_classic()
  
  
  
  ##GAM between size and extracted residuals
  ##extract residuals
  EBSf_res_orig<- residuals(EBSf_size_gam)
  #absolute values only
  EBSf_res_size <- EBSf_res_orig^2
  EBSf <- EBSf %>% mutate(res_oto_size=EBSf_res_size)
  
  #The GAM of size against the extracted residuals
  EBSf_res_gam <- gam(EBSf_res_size~s(EBSf$prev))
  gam.check(EBSf_res_gam)
  summary(EBSf_res_gam)
  
  EBSf <- EBSf %>% mutate(res_preds=predict(EBSf_res_gam))
  ggplot(EBSf,aes(x=prev, y=res_oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("EBSf")+
    theme_classic()
  
  ##now for ETAS males
  ETASm_size_gam <- gam(ETASm$oto_size~s(ETASm$prev, k=4))
  ##using gam.check helps to influence the choice of the value of k
  ##4 seems to be the best to use in this case
  gam.check(ETASm_size_gam)
  
  summary(ETASm_size_gam)
  
  #R-sq=0.978
  #Deviance explained = 97.8%
  #GCV=0.0018687
  #GCV is the measure of the degree of smoothness of the function
  ETASm <- ETASm %>% mutate(preds=predict(ETASm_size_gam))
  
  ggplot(ETASm,aes(x=prev, y=oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
    ylab("s'")+
    xlab("s")+
    ggtitle("ETASm")+
    theme_classic()
  
  
  ##GAM between size and extracted residuals
  ##extract residuals
  ETASm_res_orig<- residuals(ETASm_size_gam)
  #absolute values only
  ETASm_res_size <- (ETASm_res_orig)^2
  ETASm <- ETASm %>% mutate(res_oto_size=ETASm_res_size)
  
  #The GAM of size against the extracted residuals
  ETASm_res_gam <- gam(ETASm_res_size~s(ETASm$prev))
  gam.check(ETASm_res_gam)
  summary(ETASm_res_gam)
  
  ETASm <- ETASm %>% mutate(res_preds=predict(ETASm_res_gam))
  ggplot(ETASm,aes(x=prev, y=res_oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("ETASm")+
    theme_classic()
  
  
  ##now for ETAS females
  ETASf_size_gam <- gam(ETASf$oto_size~s(ETASf$prev, k=4))
  ##using gam.check helps to influence the choice of the value of k
  ##4 seems to be the best to use in this case
  gam.check(ETASf_size_gam)
  
  summary(ETASf_size_gam)
  
  #R-sq=0.978
  #Deviance explained = 97.8%
  #GCV=0.0018687
  #GCV is the measure of the degree of smoothness of the function
  ETASf <- ETASf %>% mutate(preds=predict(ETASf_size_gam))
  
  ggplot(ETASf,aes(x=prev, y=oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
    ylab("s'")+
    xlab("s")+
    ggtitle("ETASf")+
    theme_classic()
  
  
  
  ##GAM between size and extracted residuals
  ##extract residuals
  ETASf_res_orig<- residuals(ETASf_size_gam)
  #absolute values only
  ETASf_res_size <- (ETASf_res_orig)^2
  ETASf <- ETASf %>% mutate(res_oto_size=ETASf_res_size)
  
  #The GAM of size against the extracted residuals
  ETASf_res_gam <- gam(ETASf_res_size~s(ETASf$prev))
  gam.check(ETASf_res_gam)
  summary(ETASf_res_gam)
  
  ETASf <- ETASf %>% mutate(res_preds=predict(ETASf_res_gam))
  ggplot(ETASf,aes(x=prev, y=res_oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("ETASf")+
    theme_classic()
  
  ##now for NSW males
  NSWm_size_gam <- gam(NSWm$oto_size~s(NSWm$prev, k=4))
  ##using gam.check helps to influence the choice of the value of k
  ##4 seems to be the best to use in this case
  gam.check(NSWm_size_gam)
  
  summary(NSWm_size_gam)
  
  #R-sq=0.978
  #Deviance explained = 97.8%
  #GCV=0.0018687
  #GCV is the measure of the degree of smoothness of the function
  NSWm <- NSWm %>% mutate(preds=predict(NSWm_size_gam))
  ggplot(NSWm,aes(x=prev, y=oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
    ylab("s'")+
    xlab("s")+
    ggtitle("NSWm")+
    theme_classic()
  
  
  ##GAM between size and extracted residuals
  ##extract residuals
  NSWm_res_orig<- residuals(NSWm_size_gam)
  #absolute values only
  NSWm_res_size <-(NSWm_res_orig)^2
  NSWm <- NSWm %>% mutate(res_oto_size=NSWm_res_size)
  
  #The GAM of size against the extracted residuals
  NSWm_res_gam <- gam(NSWm_res_size~s(NSWm$prev))
  gam.check(NSWm_res_gam)
  summary(NSWm_res_gam)
  
  NSWm <- NSWm %>% mutate(res_preds=predict(NSWm_res_gam))
  ggplot(NSWm,aes(x=prev, y=res_oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("NSWm")+
    theme_classic()
  
  
  ##now for NSW females
  NSWf_size_gam <- gam(NSWf$oto_size~s(NSWf$prev, k=4))
  ##using gam.check helps to influence the choice of the value of k
  ##4 seems to be the best to use in this case
  gam.check(NSWf_size_gam)
  
  summary(NSWf_size_gam)
  
  NSWf <- NSWf %>% mutate(preds=predict(NSWf_size_gam))
  ggplot(NSWf,aes(x=prev, y=oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
    ylab("s'")+
    xlab("s")+
    ggtitle("NSWf")+
    theme_classic()
  
  
  ##GAM between size and extracted residuals
  ##extract residuals
  NSWf_res_orig<- residuals(NSWf_size_gam)
  #absolute values only
  NSWf_res_size <- (NSWf_res_orig)^2
  NSWf <- NSWf %>% mutate(res_oto_size=NSWf_res_size)
  
  #The GAM of size against the extracted residuals
  NSWf_res_gam <- gam(NSWf_res_size~s(NSWf$prev))
  gam.check(NSWf_res_gam)
  summary(NSWf_res_gam)
  
  NSWf <- NSWf %>% mutate(res_preds=predict(NSWf_res_gam))
  ggplot(NSWf,aes(x=prev, y=res_oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("NSWf")+
    theme_classic()
  
  ##now for WTAS males
  WTASm_size_gam <- gam(WTASm$oto_size~s(WTASm$prev, k=4))
  ##using gam.check helps to influence the choice of the value of k
  ##4 seems to be the best to use in this case
  gam.check(WTASm_size_gam)
  
  summary(WTASm_size_gam)
  
  #R-sq=0.978
  #Deviance explained = 97.8%
  #GCV=0.0018687
  #GCV is the measure of the degree of smoothness of the function
  WTASm <- WTASm %>% mutate(preds=predict(WTASm_size_gam))
  ggplot(WTASm,aes(x=prev, y=oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
    ylab("s'")+
    xlab("s")+
    ggtitle("WTASm")+
    theme_classic()
  
  
  ##GAM between size and extracted residuals
  ##extract residuals
  WTASm_res_orig<- residuals(WTASm_size_gam)
  #absolute values only
  WTASm_res_size <- (WTASm_res_orig)^2
  WTASm <- WTASm %>% mutate(res_oto_size=WTASm_res_size)
  
  #The GAM of size against the extracted residuals
  WTASm_res_gam <- gam(WTASm_res_size~s(WTASm$prev))
  gam.check(WTASm_res_gam)
  summary(WTASm_res_gam)
  
  WTASm <- WTASm %>% mutate(res_preds=predict(WTASm_res_gam))
  ggplot(WTASm,aes(x=prev, y=res_oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("WTASm")+
    theme_classic()
  
  
  ##now for WTAS females
  WTASf_size_gam <- gam(WTASf$oto_size~s(WTASf$prev, k=4))
  ##using gam.check helps to influence the choice of the value of k
  ##4 seems to be the best to use in this case
  gam.check(WTASf_size_gam)
  
  summary(WTASf_size_gam)
  
  #R-sq=0.978
  #Deviance explained = 97.8%
  #GCV=0.0018687
  #GCV is the measure of the degree of smoothness of the function
  
  WTASf <- WTASf %>% mutate(preds=predict(WTASf_size_gam))
  ggplot(WTASf,aes(x=prev, y=oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=preds), size=1.3, colour="steelblue")+
    ylab("s'")+
    xlab("s")+
    ggtitle("WTASf")+
    theme_classic()
  
  
  ##GAM between size and extracted residuals
  ##extract residuals
  WTASf_res_orig<- residuals(WTASf_size_gam)
  #absolute values only
  WTASf_res_size <- (WTASf_res_orig)^2
  WTASf <- WTASf %>% mutate(res_oto_size=WTASf_res_size)
  
  #The GAM of size against the extracted residuals
  WTASf_res_gam <- gam(WTASf_res_size~s(WTASf$prev))
  gam.check(WTASf_res_gam)
  summary(WTASf_res_gam)
  
  WTASf <- WTASf %>% mutate(res_preds=predict(WTASf_res_gam))
  ggplot(WTASf,aes(x=prev, y=res_oto_size))+
    geom_point(size=1)+
    geom_line(aes(x=prev, y=res_preds), size=1.3, colour="steelblue")+
    ylab("r^2")+
    xlab("s")+
    ggtitle("WTASf")+
    theme_classic()
  
  
 
  