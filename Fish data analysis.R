#Fish data analysis
#11/12/2018

#Setting up-----------
#Libraries required
library(dplyr)
library(ggplot2)

#Data required
fish <- read.csv("./otoliths (working)/data_derived/data_otolith_complete.csv")

#See the Data sleuthing notes Rmarkdown file for some general comments about 
#the datasets and its constituent variables

#explore the data a bit
glimpse(fish)

str(fish)

#make a few plots of variables of interest- start with the effect of temperature on some variables
ggplot(fish, aes(x=bottomtemp1, y=Increment))+
  geom_point(size=1)

ggplot(fish, aes(x=bottomtemp1, y=log_incr))+
  geom_point(size=1)

ggplot(fish, aes(x=bottomtemp1, y=growth))+
  geom_point(size=1)

#Shows that as the temperature increases so does the 
#population-wide average growth

#for radius and fish length must group by fish ID
ggplot(fish, aes(x=bottomtemp1, y=radius))+
  geom_point(size=1)


#make a few plots of variables of interest- have a look at age in years now

ggplot(fish, aes(x=Age, y=growth))+
  geom_point(size=1)

ggplot(fish, aes(x=Age, y=radius))+
  geom_point(size=1)

ggplot(fish, aes(x=Age, y=Increment))+
  geom_point(size=1)

ggplot(fish, aes(x=Age, y=log_incr))+
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

#Effect of area caught

ggplot(fish, aes(x=area, y=growth))+
  geom_boxplot()

ggplot(fish, aes(x=area, y=radius))+
  geom_boxplot()

ggplot(fish, aes(x=area, y=Increment))+
  geom_boxplot()

ggplot(fish, aes(x=area, y=log_incr))+
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

#plotting fish length against otolith size variables

ggplot(fish, aes(x=floorlength, y=growth))+
  geom_point(size=1)

ggplot(fish, aes(x=floorlength, y=radius))+
  geom_point(size=1)

ggplot(fish, aes(x=floorlength, y=Increment))+
  geom_point(size=1)

ggplot(fish, aes(x=floorlength, y=log_incr))+
  geom_point(size=1)

#Fishing gear effects

ggplot(fish, aes(x=gear, y=growth))+
  geom_boxplot()

ggplot(fish, aes(x=gear, y=radius))+
  geom_boxplot()

ggplot(fish, aes(x=gear, y=Increment))+
  geom_boxplot()

ggplot(fish, aes(x=gear, y=log_incr))+
  geom_boxplot()
