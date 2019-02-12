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
not_1 <- mutate(not_1, prev=fish[ind_not_1-1,34])
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
EBSm_size_gam <- gam(EBSm$oto_size~s(EBSm$prev, k=4))
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
EBSm_res_size<- residuals(EBSm_size_gam)
#absolute values only
EBSm_res_size <- (EBSm_res_size)^2
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
EBSf_res_size<- residuals(EBSf_size_gam)
#absolute values only
EBSf_res_size <- (EBSf_res_size)^2
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
ETASm_res_size<- residuals(ETASm_size_gam)
#absolute values only
ETASm_res_size <- (ETASm_res_size)^2
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
ETASf_res_size<- residuals(ETASf_size_gam)
#absolute values only
ETASf_res_size <- (ETASf_res_size)^2
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
NSWm_res_size<- residuals(NSWm_size_gam)
#absolute values only
NSWm_res_size <-(NSWm_res_size)^2
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
NSWf_res_size<- residuals(NSWf_size_gam)
#absolute values only
NSWf_res_size <- (NSWf_res_size)^2
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
WTASm_res_size<- residuals(WTASm_size_gam)
#absolute values only
WTASm_res_size <- (WTASm_res_size)^2
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
WTASf_res_size<- residuals(WTASf_size_gam)
#absolute values only
WTASf_res_size <- (WTASf_res_size)^2
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









##Now we can move on to the transformations of these
#This shows that as size increases the residuals are decreasing 
#We would like to remove this effect as much as possible
#So we will try box-cox transformations

##Box-cox transformations of the otolith size data

##I will try the transformations out on the NSWm data first as it has a particularly noticeable 
##relationship

#Transform the otolith size data using bcPower from the 'car' package
NSWm <- mutate(NSWm, NSWm_bc_oto_size=bcPower(oto_size, lambda = seq(-2, 2, length=50)))
NSWm <- mutate(NSWm, NSWm_bc_prev=bcPower(prev, lambda = seq(-2,2, length=50)))

#create the box-cox GAM
NSWm_bc_size_gam <- gam(NSWm$NSWm_bc_oto_size~s(NSWm$NSWm_bc_prev, k=4))
gam.check(NSWm_bc_size_gam)
summary(NSWm_bc_size_gam)
NSWm <- NSWm %>% mutate(bc_preds=predict(NSWm_bc_size_gam))
ggplot(NSWm,aes(x=NSWm_bc_prev, y=NSWm_bc_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=NSWm_bc_prev, y=bc_preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("Box-Cox NSWm")+
  theme_classic()

#Extract the residuals from this
NSWm_res_bc_size <- residuals(NSWm_bc_size_gam)
NSWm_res_bc_size <- (NSWm_res_bc_size)^2
NSWm <- NSWm %>% mutate(res_bc_size=NSWm_res_bc_size)

#GAM of the residuals
NSWm_res_bc_gam <- gam(NSWm$res_bc_size~s(NSWm$NSWm_bc_prev))
gam.check(NSWm_res_bc_gam)
summary(NSWm_res_bc_gam)

NSWm <- NSWm %>% mutate(res_bc_preds=predict(NSWm_res_bc_gam))
ggplot(NSWm,aes(x=NSWm_bc_prev, y=res_bc_size))+
  geom_point(size=1)+
  geom_line(aes(x=NSWm_bc_prev, y=res_bc_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Box-Cox NSWm")+
  theme_classic()

##That didn't seem to work
##This method seems to create a lot of NAs within the dataset

#same relationship found- try another box-cox rather than bcPower (try 'powerTransform' and 'boxCoxVariable'- all from 'car' package)
#Transform the otolith size data
NSWm <- mutate(NSWm, bcv_oto_size=boxCoxVariable(oto_size))
NSWm <- mutate(NSWm, bcv_oto_prev=boxCoxVariable(prev))

#create the boxcoxVariable GAM extract residuals and then plot the residual GAM
NSWm_bcv_size_gam <- gam(NSWm$bcv_oto_size~s(NSWm$bcv_oto_prev, k=4))
gam.check(NSWm_bcv_size_gam)
summary(NSWm_bcv_size_gam)
NSWm <- NSWm %>% mutate(bcv_preds=predict(NSWm_bcv_size_gam))
ggplot(NSWm,aes(x=bcv_oto_prev, y=bcv_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=bcv_oto_prev, y=bcv_preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("Box-Cox Variable NSWm")+
  theme_classic()

#Extract the residuals from this
NSWm_res_bcv_size <- residuals(NSWm_bcv_size_gam)
NSWm_res_bcv_size <- (NSWm_res_bcv_size)^2
NSWm <- NSWm %>% mutate(res_bcv_size=NSWm_res_bcv_size)

#GAM of the residuals
NSWm_res_bcv_gam <- gam(NSWm$res_bcv_size~s(NSWm$bcv_oto_prev))
gam.check(NSWm_res_bcv_gam)
summary(NSWm_res_bcv_gam)

NSWm <- NSWm %>% mutate(res_bcv_preds=predict(NSWm_res_bcv_gam))
ggplot(NSWm,aes(x=bcv_oto_prev, y=res_bcv_size))+
  geom_point(size=1)+
  geom_line(aes(x=bcv_oto_prev, y=res_bcv_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Box-Cox Variable NSWm")+
  theme_classic()

#Better than bcPower function-but still not right


##Try a different method for Box-Cox transformation
#The manual method 
Box=boxcox(NSWm$oto_size~1,
           lambda= seq(-6, 6, 0.1)
           )

Cox=data.frame(Box$x, Box$y)
Cox2=Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]

#To obtain the appropriate value of lambda
lambda = Cox2[1, "Box.x"]
oto_box = (NSWm$oto_size ^ lambda -1)/lambda
oto_box_prev = (NSWm$prev^lambda-1)/lambda

#Make a GAM of this new boxcox data


#create the box-cox GAM
NSWm_box_size_gam <- gam(oto_box~s(oto_box_prev, k=4))
gam.check(NSWm_box_size_gam)
summary(NSWm_box_size_gam)
NSWm <- NSWm %>% mutate(box_preds=predict(NSWm_box_size_gam))
ggplot(NSWm,aes(x=oto_box_prev, y=oto_box))+
  geom_point(size=1)+
  geom_line(aes(x=oto_box_prev, y=box_preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("Manual Box-Cox NSWm")+
  theme_classic()



#Extract the residuals from this
NSWm_res_box_size <- residuals(NSWm_box_size_gam)
NSWm_res_box_size <- (NSWm_res_box_size)^2
NSWm <- NSWm %>% mutate(res_box_size=NSWm_res_box_size)

#GAM of the residuals
NSWm_res_box_gam <- gam(NSWm$res_box_size~s(oto_box_prev, k=4))
gam.check(NSWm_res_box_gam)
summary(NSWm_res_box_gam)

NSWm <- NSWm %>% mutate(res_box_preds=predict(NSWm_res_box_gam))
ggplot(NSWm,aes(x=oto_box_prev, y=res_box_size))+
  geom_point(size=1)+
  geom_line(aes(x=oto_box_prev, y=res_box_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox NSWm")+
  theme_classic()

#No still not good


##Try boxCoxVariable for the ETASm dataset now to see what happens
ETASm <- mutate(ETASm, bcv_oto_size=boxCoxVariable(oto_size))
ETASm <- mutate(ETASm, bcv_oto_prev=boxCoxVariable(prev))

#create the boxcoxVariable GAM extract residuals and then plot the residual GAM
ETASm_bcv_size_gam <- gam(ETASm$bcv_oto_size~s(ETASm$bcv_oto_prev, k=4))
gam.check(ETASm_bcv_size_gam)
summary(ETASm_bcv_size_gam)
ETASm <- ETASm %>% mutate(bcv_preds=predict(ETASm_bcv_size_gam))
ggplot(ETASm,aes(x=bcv_oto_prev, y=bcv_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=bcv_oto_prev, y=bcv_preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("Box-Cox Variable ETASm")+
  theme_classic()

#Extract the residuals from this
ETASm_res_bcv_size <- residuals(ETASm_bcv_size_gam)
ETASm_res_bcv_size <- (ETASm_res_bcv_size)^2
ETASm <- ETASm %>% mutate(res_bcv_size=ETASm_res_bcv_size)

#GAM of the residuals
ETASm_res_bcv_gam <- gam(ETASm$res_bcv_size~s(ETASm$bcv_oto_prev))
gam.check(ETASm_res_bcv_gam)
summary(ETASm_res_bcv_gam)

ETASm <- ETASm %>% mutate(res_bcv_preds=predict(ETASm_res_bcv_gam))
ggplot(ETASm,aes(x=bcv_oto_prev, y=res_bcv_size))+
  geom_point(size=1)+
  geom_line(aes(x=bcv_oto_prev, y=res_bcv_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Box-Cox Variable ETASm")+
  theme_classic()

##No that is not a good one
##Try the manual method now for ETASm

ETASm_Box=boxcox(ETASm$oto_size~1,
           lambda= seq(-6, 6, 0.1)
)

ETASm_Cox=data.frame(Box$x, Box$y)
ETASm_Cox2=ETASm_Cox[with(Cox, order(-Cox$Box.y)),]
ETASm_Cox2[1,]

#To obtain the appropriate value of lambda
lambda = ETASm_Cox2[1, "Box.x"]
ETASm_oto_box = (ETASm$oto_size ^ lambda -1)/lambda
ETASm_oto_box_prev = (ETASm$prev^lambda-1)/lambda

#Make a GAM of this new boxcox data


#create the box-cox GAM
ETASm_box_size_gam <- gam(ETASm_oto_box~s(ETASm_oto_box_prev, k=4))
gam.check(ETASm_box_size_gam)
summary(ETASm_box_size_gam)
ETASm <- ETASm %>% mutate(box_preds=predict(ETASm_box_size_gam))
ggplot(ETASm,aes(x=ETASm_oto_box_prev, y=ETASm_oto_box))+
  geom_point(size=1)+
  geom_line(aes(x=ETASm_oto_box_prev, y=box_preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("Manual Box-Cox ETASm")+
  theme_classic()



#Extract the residuals from this
ETASm_res_box_size <- residuals(ETASm_box_size_gam)
ETASm_res_box_size <- (ETASm_res_box_size)^2
ETASm <- ETASm %>% mutate(res_box_size=ETASm_res_box_size)

#GAM of the residuals
ETASm_res_box_gam <- gam(ETASm$res_box_size~s(ETASm_oto_box_prev, k=4))
gam.check(ETASm_res_box_gam)
summary(ETASm_res_box_gam)

ETASm <- ETASm %>% mutate(res_box_preds=predict(ETASm_res_box_gam))
ggplot(ETASm,aes(x=ETASm_oto_box_prev, y=res_box_size))+
  geom_point(size=1)+
  geom_line(aes(x=ETASm_oto_box_prev, y=res_box_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox ETASm")+
  theme_classic()

##No not good