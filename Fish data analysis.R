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
EBSm_res_size <- abs(EBSm_res_size)^2
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
EBSf_res_size <- abs(EBSf_res_size)^2
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
ETASm_res_size <- abs(ETASm_res_size)^2
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
ETASf_res_size <- abs(ETASf_res_size)^2
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
NSWm_res_size <- abs(NSWm_res_size)^2
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
NSWf_res_size <- abs(NSWf_res_size)^2
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
WTASm_res_size <- abs(WTASm_res_size)^2
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
WTASf_res_size <- abs(WTASf_res_size)^2
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
NSWm_res_bc_size <- abs(NSWm_res_bc_size)^2
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

#same relationship found- try another box-cox rather than bcPower (try 'powerTransform' and 'boxCoxVariable' and 'boxTidwell'- all from 'car' package)
#Transform the otolith size data
fish <- mutate(fish, bcv_oto_size=boxCoxVariable(oto_size))

#create the boxcoxVariable GAM extract residuals and then plot the residual GAM
bcv_size_gam <- gam(fish$bcv_oto_size[ind_not_1]~s(fish$bcv_oto_size[ind_not_1-1], k=4, by=not_1$zone))
gam.check(bcv_size_gam)
summary(bcv_size_gam)
plot(bcv_size_gam)

res_bcv_size <- residuals(bcv_size_gam)
res_bcv_size <- abs(res_bcv_size)

res_bcv_gam <- gam(res_bcv_size~s(fish$bcv_oto_size[ind_not_1-1], by=not_1$zone))
gam.check(res_bcv_gam)
summary(res_bcv_gam)
plot(res_bcv_gam)
#Different looking relationship from the bcPower function-but still not right- not too bad compared with the other one
#The one for ETAS looks decent


#I'll try a log transformation now to see what happens for the GAMs

log_size_gam <- gam(fish$log_oto_size[ind_not_1]~s(fish$log_oto_size[ind_not_1-1], k=4, by=not_1$zone))
gam.check(log_size_gam)
summary(log_size_gam)
plot(log_size_gam)

res_log_size <- residuals(log_size_gam)
res_log_size <- abs(res_log_size)

res_log_gam <- gam(res_log_size~s(fish$log_oto_size[ind_not_1-1], by=not_1$zone))
gam.check(res_log_gam)
summary(res_log_gam)
plot(res_log_gam)
##no still not good- similar relationship as seen with bcPower function

##Try a different method for Box-Cox transformation
Box=boxcox(NSWm$oto_size~1,
           lambda= seq(-6, 6, 0.1)
           )

Cox=data.frame(Box$x, Box$y)
Cox2=Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]

#To obtain the appropriate value of lambda
lambda = Cox2[1, "Box.x"]
oto_box = (oto_size ^ lambda -1)/lambda

#Make a GAM of this new boxcox data

box_size_gam <- gam(oto_box[ind_not_1]~s(oto_box[ind_not_1-1], k=4, by=not_1$zone))
gam.check(box_size_gam)
summary(box_size_gam)
plot(box_size_gam)

res_box_size <- residuals(box_size_gam)
res_box_size <- abs(res_box_size)

res_box_gam <- gam(res_box_size~s(oto_box[ind_not_1-1], by=not_1$zone))
gam.check(res_box_gam)
summary(res_box_gam)
plot(res_box_gam)
#No still not good
