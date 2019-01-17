##Making Fede's code work

library(mgcv)
library(dplyr)
library(ggplot2)
library(MASS)
library(plot3D)
library(MuMIn)

######## Data arrangement
dd<-read.csv("./otoliths (working)/data_derived/data_otolith_complete.csv")
##Making a new coloumn called otolith size and then logging that and making a
##new csv from that- cqn't get this to work
oto_size<-rep(0,27429) 
for (i in 1:27429) oto_size[i] <- sum(Increment[((i-Age[i])+1):i])
dd<-data.frame(cbind(dd,oto_size))
dd<-mutate(dd, log_oto_size=log(oto_size))
write.csv(dd,"extended_data_oto.csv")
########

####### Parallel data oto
dd<-read.csv("extended_data_oto.csv")
attach(dd)
##which fish have an age that isn't 1

ind_not_1<-which(Age!=1)
new<-filter(dd, Age!=1)
##a new column added called prev which consists of the age not 1 data
#Not sure what the 1 and 33 do and why these numbers are chosen
new<-mutate(new, prev=dd[ind_not_1-1,33])
write.csv(new,"parallel_data_oto.csv")
#######

dd<-read.csv("extended_data_oto.csv")
#again don't know what is happening here
dd<-dd[, 4:34]

# At the very beginning I conduct the analysis only on a region.
dp<-dd[which(dd[,1]=="NSW"),]
attach(dp)

# Let's start to look into the size-dependant growth 
ind_not_1<-which(Age!=1)


##############################################################################
### Preliminary Analysis: normal  and generalized linear regression models ###
##############################################################################

# log_Increment(t+1) vs. log_Size (t)     0.5752
windows()
par(mfrow=c(2,2))
plot(log_oto_size[ind_not_1-1], log_incr[ind_not_1], pch=3,
     xlab="log_Size(t)", ylab="log_Increment(t+1)")
fit<-lm(log_incr[ind_not_1]~log_oto_size[ind_not_1-1])
abline(fit$coef, col='red', lwd=2)
plot(fit,which=1:3)
summary(fit)

# Increment(t+1) vs. Size (t)     0.5445 no buono
windows()
par(mfrow=c(2,2))
plot(oto_size[ind_not_1-1], Increment[ind_not_1], pch=3,
     xlab="Size(t)", ylab="Increment(t+1)")
fit<-lm(Increment[ind_not_1]~oto_size[ind_not_1-1])
abline(fit$coef, col='red', lwd=2)
plot(fit,which=1:3)
summary(fit)


# Increment(t+1) vs. Increment (t)  Disgustoso
windows()
par(mfrow=c(2,2))
plot(Increment[ind_not_1-1], Increment[ind_not_1], pch=3,
     xlab="Increment(t)", ylab="Increment(t+1)")
fit<-lm(Increment[ind_not_1]~Increment[ind_not_1-1])
abline(fit$coef, col='red', lwd=2)
plot(fit,which=1:3)
summary(fit)

# log_incr(t+1) vs. log_incr(t)  0.6151
windows()
par(mfrow=c(2,2))
plot(log_incr[ind_not_1-1], log_incr[ind_not_1], pch=3,
     xlab="log_incr(t)", ylab="log_incr(t+1)")
fit<-lm(log_incr[ind_not_1]~log_incr[ind_not_1-1])
abline(fit$coef, col='red', lwd=2)
plot(fit,which=1:3)
summary(fit)
# Add age to this one, that is the best model. 
windows()
par(mfrow=c(2,2))
plot(log_incr[ind_not_1-1], log_incr[ind_not_1], pch=3,
     xlab="log_incr(t)", ylab="log_incr(t+1)")
fit1<-lm(log_incr[ind_not_1]~log_incr[ind_not_1-1]+Age[ind_not_1])
abline(fit$coef, col='red', lwd=2)
lines(sort(log_incr[ind_not_1-1]), fit1$fitted.values[order(log_incr[ind_not_1-1])], 
      col='blue', lwd=2)
data.prov<-sort(log_incr[ind_not_1-1])
plot(fit1,which=1:3)
summary(fit1)
# R-squared:  0.6495. Diagnostics plot are worse than only R. That could be because of the
# collinearity between the two regressors. 


#### GLM
# What if trying gamma distribution? No good results on logarithm
fit.glm<-glm(abs(log_incr[ind_not_1])~abs(log_incr[ind_not_1-1])+Age[ind_not_1], family="Gamma")
summary(fit.glm)
plot(fit.glm)
r.squaredLR(fit.glm) # R2=0.59
windows()
par(mfrow=c(2,2))
plot(abs(log_incr[ind_not_1-1]), abs(log_incr[ind_not_1]), pch=3,
     xlab="abs log_incr(t)", ylab="abs log_incr(t+1)")
lines(sort(abs(log_incr[ind_not_1-1])), fit.glm$fitted.values[order(abs(log_incr[ind_not_1-1]))], 
      col='blue', lwd=2)
plot(fit.glm,which=1:3)
# Not good diagnostics (R2=0.59)

# GLM on increment
fit.glm1<-glm(Increment[ind_not_1]~Increment[ind_not_1-1]+Age[ind_not_1], family="Gamma")
summary(fit.glm1)
r.squaredLR(fit.glm1) # R2=0.62
windows()
par(mfrow=c(2,2))
plot(Increment[ind_not_1-1], Increment[ind_not_1], pch=3,
     xlab="Increment(t)", ylab="Increment(t+1)")
lines(sort(Increment[ind_not_1-1]), fit.glm1$fitted.values[order(Increment[ind_not_1-1])], 
      col='blue', lwd=2)
plot(fit.glm1,which=1:3)
# This seems to be the best model, along with the LM on logarithm. Diagnostics plots are
# good and the amount of variability explained is big. 


############
#### GAM ###
############

### Consider now Generalized Additive Models (Size(t+1) vs Size(t) or its logarithm)
?gam

# Reference: linear model
fit<-lm(oto_size[ind_not_1]~oto_size[ind_not_1-1])
summary(fit)

# Ok let's start. 
fit.gam<-gam(oto_size[ind_not_1]~s(oto_size[ind_not_1-1], k=4))
gam.check(fit.gam)
summary(fit.gam) # R-sq.(adj) =  0.952, GCV = 0.0025867
# Automatically R choose rank=10
# I tried to manual select the knots position, no improvements.
plot(fit.gam)

plot(oto_size[ind_not_1-1],oto_size[ind_not_1], pch=3, xlab="Otolith size (t)",
     ylab="Otolith size (t+1)")
abline(0,1, col='red', lwd=2)
lines(sort(oto_size[ind_not_1-1]), fit.gam$fitted.values[order(oto_size[ind_not_1-1])], 
      col='blue', lwd=2)
# The result is good: the estimated growth curve is asymptotic to y=x+a, because of the
# horizontal asymptote in otolith increment width

# Try to insert the temperature factor in the model
fit.gam1<-gam(oto_size[ind_not_1]~s(oto_size[ind_not_1-1])+bottomtemp1[ind_not_1])
gam.check(fit.gam1)
summary(fit.gam1)
# The temperature does not produce a better explanation (R-sq.(adj) =  0.953)
# It has an effect, but the predictive capacity does not increase a lot. GCV decreases
# from 0.0025867 to 0.0025187 (3%)

# Temperature as a smooth effect
fit.gam11<-gam(oto_size[ind_not_1]~s(oto_size[ind_not_1-1])+s(bottomtemp1[ind_not_1]))
gam.check(fit.gam11)
summary(fit.gam11)
# GCV 0.0024923
plot(fit.gam11)
# Temperature seems not to influence otolith growth. 

# Incr(t) + AGe effect  THIS IS THE MODEL!! 
fit.gam2<-gam(oto_size[which(Age!=1)]~s(oto_size[which(Age!=1)-1])+Age[which(Age!=1)])
gam.check(fit.gam2)
summary(fit.gam2)
plot(fit.gam2)

plot(oto_size[ind_not_1-1],oto_size[ind_not_1], pch=3, xlab="Otolith size (t)",
     ylab="Otolith size (t+1)")
abline(0,1, col='red', lwd=2)
lines(sort(oto_size[ind_not_1-1]), fit.gam2$fitted.values[order(oto_size[ind_not_1-1])], 
      col='blue', lwd=2)
# Ok better! GCV = 0.002371, R-sq.(adj) =  0.957 and the plot (if Age considered smooth) is good

# increment(t) + Smooth AGe effect - No substantial difference
fit.gam22<-gam(oto_size[which(Age!=1)]~s(oto_size[which(Age!=1)-1])+s(Age[which(Age!=1)]))
summary(fit.gam22)
plot(fit.gam22)
plot(oto_size[ind_not_1-1],oto_size[ind_not_1], pch=3, xlab="Otolith size (t)",
     ylab="Otolith size (t+1)")
abline(0,1, col='red', lwd=2)
lines(sort(oto_size[ind_not_1-1]), fit.gam22$fitted.values[order(oto_size[ind_not_1-1])], 
      col='blue', lwd=2)

# Only Age
fit.gam3<-gam(oto_size[ind_not_1]~s(Age[ind_not_1]))
summary(fit.gam3)
plot(fit.gam3)
# Interesting: if Age is considered alone its R-sq.(adj) =  0.794. In fact it exists
# a correlation between Age and Size, the two regressors are then collinear. 

# I did a "forced autoregressive model", then the predictive function should be builded. 


### LOGARITHMIC RESPONSE - no better
#####
fit.gaml<-gam(log_oto_size[ind_not_1]~s(log_oto_size[ind_not_1-1]))
gam.check(fit.gaml)
summary(fit.gaml)
# Automatically R choose rank=10
# I tried to manual select the knots position, no improvements.
plot(fit.gaml)
# GCV = 0.0037423

# With age?
fit.gaml1<-gam(log_oto_size[ind_not_1]~s(log_oto_size[ind_not_1-1])+s(Age[ind_not_1]))
plot(fit.gaml1)
# The age effect is more clear than temperature
fit.gaml2<-gam(log_oto_size[ind_not_1]~s(log_oto_size[ind_not_1-1])+Age[ind_not_1])
gam.check(fit.gaml2)
summary(fit.gaml2)
# GCV = 0.0035763
######






#########################
#### Local regression ### # Here I am wornking on the "parallel" dataset
#########################
dd<-read.csv("parallel_data_oto.csv")
dd<-data[,5:36]
data<-dd[which(zone=="NSW"),]
attach(data)

# Tricubic kernel 
fit.loess<-loess(oto_size~prev, span = .008)
summary(fit.loess)
plot(prev,oto_size, pch=3, xlab="Otolith size (t)",
     ylab="Otolith size (t+1)")
abline(0,1, col='red', lwd=2)
lines(sort(prev), fit.loess$fitted[order(prev)],
      type='l',col='blue', lwd=2)
SStot<-sum((oto_size-mean(oto_size))^2)
SSres<-sum((fit.loess$residuals)^2)
R2<-1-SSres/SStot
R2
# We obtain a 96.2 r-squared, with the best asymptotic behaviour 
# (asymptotic to y=x+a)
# Add Age
fit.loess2<-loess(oto_size[ind_not_1]~oto_size[ind_not_1-1]+Age[ind_not_1], span = 0.01)
summary(fit.loess2)
plot(oto_size[ind_not_1-1],oto_size[ind_not_1], pch=3, xlab="Otolith size (t)",
     ylab="Otolith size (t+1)")
abline(0,1, col='red', lwd=2)
lines(sort(oto_size[ind_not_1-1]), sort(fit.loess2$fitted),
      type='l',col='blue', lwd=2)
SStot<-sum((oto_size[ind_not_1]-mean(oto_size[ind_not_1]))^2)
SSres<-sum((fit.loess2$residuals)^2)
R2<-1-SSres/SStot
R2
# R2=96.6, the same as GAM

# Nearest neighbours
res<-rep(NA,100)
for (k in 20:40) {   # Consider the k-nearest neighbours in fitting the model 
  k
  bw<-k
  fitted<-NULL
  real<-NULL
  for (i in 1:dim(data)[1]) {
    sapply(1:dim(data)[1], FUN = function(i){
      finestra<-data[order(abs(prev - prev[i]))[1:k],] 
      fit<-lm(oto_size~prev+Age, data=finestra[-1,])
      fitted<-c(fitted,predict.lm(fit, finestra[1,]))
      real<-c(real,finestra$oto_size[1])
      return(sum((fitted-real)^2)/length(fitted))
    })
    res[k]<-sum((fitted-real)^2)/length(fitted)
  }
  plot(20:40,res[20:40])
  
  # Gaussian Kernel
  
  # Epanechnikov
  
  
  ####################################
  #### Svalvolata on the road!!! #####
  ####################################
  
  den<-kde2d(oto_size[ind_not_1-1],oto_size[ind_not_1], n=50)
  plot(oto_size[ind_not_1-1],oto_size[ind_not_1], pch=3, xlab="Otolith size (t)",
       ylab="Otolith size (t+1)")
  abline(0,1,col='red',lwd=2)
  contour(den, drawlabels=FALSE, nlevels=8, add=TRUE, col='darkorange',lwd=1.5) 
  persp3D(z=den$z, phi = 45, theta = 45)
  
  denl<-kde2d(log_oto_size[ind_not_1-1],log_oto_size[ind_not_1], n=50)
  plot(log_oto_size[ind_not_1-1],log_oto_size[ind_not_1], pch=3, xlab="Log Otolith size (t)",
       ylab="Log Otolith size (t+1)")
  abline(0,1,col='red',lwd=2)
  contour(denl, drawlabels=FALSE, nlevels=8, add=TRUE, col='darkorange',lwd=1.5)
  persp3D(z=denl$z, phi = 45, theta = 45)