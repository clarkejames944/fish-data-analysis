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
not_1 <- mutate(not_1, prev=fish[ind_not_1-1,33])
write.csv(not_1,"parallel_data_oto.csv")

##Now all is set up to do a GAM
##Oto_size against previous year
?gam

size_gam <- gam(oto_size[ind_not_1]~s(oto_size[ind_not_1-1], k=4, by=not_1$zone))
##using gam.check helps to influence the choice of the value of k
##4 seems to be the best to use in this case
gam.check(size_gam)

summary(size_gam)
#Number of observations was 19130
#This means that the analysis doesn't include the age 1 fish 
#R-sq=0.978
#Deviance explained = 97.8%
#GCV=0.0018687
#GCV is the measure of the degree of smoothness of the function

plot(size_gam)


##GAM between size and extracted residuals
##extract residuals
res_size<- residuals(size_gam)
#absolute values only
res_size <- abs(res_size)

#The GAM of size against the extracted residuals
res_gam <- gam(res_size~s(oto_size[ind_not_1-1], by=not_1$zone))
gam.check(res_gam)
summary(res_gam)

plot(res_gam)

#This shows that as size increases the residuals are decreasing 
#We would like to remove this effect as much as possible
#So we will try box-cox transformations

##Box-cox transformations of the otolith size data

#Transform the otolith size data using bcPower from the 'car' package
fish <- mutate(fish, bc_oto_size=bcPower(oto_size, lambda=TRUE))

#create the box-cox GAM
bc_size_gam <- gam(fish$bc_oto_size[ind_not_1]~s(fish$bc_oto_size[ind_not_1-1], k=4, by=not_1$zone))
gam.check(bc_size_gam)
summary(bc_size_gam)
plot(bc_size_gam)

#Extract the residuals from this
res_bc_size <- residuals(bc_size_gam)
res_bc_size <- abs(res_bc_size)

#GAM of the residuals
res_bc_gam <- gam(res_bc_size~s(fish$bc_oto_size[ind_not_1-1], by=not_1$zone))
gam.check(res_bc_gam)
summary(res_bc_gam)
plot(res_bc_gam)

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
Box=boxcox(oto_size~1,
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
