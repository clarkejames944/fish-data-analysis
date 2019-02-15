##I will try the transformations out on the NSWm data first as it has a particularly noticeable 
##relationship
##boxcox function from MASS package can only be used for lm and aov. So not that one.
#powerTransform from car package could be used
?powerTransform

#try bcPower
?bcPower

#Transform the otolith size data using bcPower from the 'car' package
NSWm <- mutate(NSWm, NSWm_bc_oto_size=bcPower(oto_size, lambda = seq(-2, 2, length=40)))
NSWm <- mutate(NSWm, NSWm_bc_prev=bcPower(prev, lambda = seq(-2,2, length=40)))

#create the box-cox GAM
NSWm_bc_size_gam <- gam(NSWm$NSWm_bc_oto_size~s(NSWm$NSWm_bc_prev, k=4))
gam.check(NSWm_bc_size_gam)
summary(NSWm_bc_size_gam)
NSWm <- NSWm %>% mutate(bc_preds=predict(NSWm_bc_size_gam))
ggplot(NSWm,aes(x=NSWm_bc_prev, y=NSWm_bc_oto_size))+
  geom_point(size=1)+
  geom_line(aes(x=NSWm_bc_prev, y=bc_preds), size=1.3, colour="steelblue")+
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Box-Cox NSWm (bcPower)")+
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
  ggtitle("Box-Cox NSWm (bcPower)")+
  theme_classic()

##That didn't seem to work
##This method seems to create a lot of NAs within the dataset

#same relationship found- try another box-cox rather than bcPower (try'boxCoxVariable'- all from 'car' package)
?boxCoxVariable
##Says that it should only be used for the response variable

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
  ylab("Box Cox s'")+
  xlab("Box Cox s")+
  ggtitle("Box-Cox Variable NSWm (boxCoxVariable)")+
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
  ylab("Box Cox r^2")+
  xlab("Box Cox s")+
  ggtitle("Box-Cox Variable NSWm (boxCoxVariable)")+
  theme_classic()

#Better than bcPower function-but still not right

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

ETASm_Cox=data.frame(ETASm_Box$x, ETASm_Box$y)
ETASm_Cox2=ETASm_Cox[with(ETASm_Cox, order(-ETASm_Cox$ETASm_Box.y)),]
ETASm_Cox2[1,]

#To obtain the appropriate value of lambda
lambda = ETASm_Cox2[1, "ETASm_Box.x"]
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

##No, still not very good

##Try the manual method now for ETASf

ETASf_Box=boxcox(ETASf$oto_size~1,
                 lambda= seq(-6, 6, 0.1)
)

ETASf_Cox=data.frame(ETASf_Box$x, ETASf_Box$y)
ETASf_Cox2=ETASf_Cox[with(ETASf_Cox, order(-ETASf_Cox$ETASf_Box.y)),]
ETASf_Cox2[1,]

#To obtain the appropriate value of lambda
lambda = ETASf_Cox2[1, "ETASf_Box.x"]
ETASf_oto_box = (ETASf$oto_size ^ lambda -1)/lambda
ETASf_oto_box_prev = (ETASf$prev^lambda-1)/lambda

#Make a GAM of this new boxcox data

#create the box-cox GAM
ETASf_box_size_gam <- gam(ETASf_oto_box~s(ETASf$prev, k=4))
gam.check(ETASf_box_size_gam)
summary(ETASf_box_size_gam)
ETASf<- ETASf %>% mutate(box_preds=predict(ETASf_box_size_gam))
ggplot(ETASf,aes(x=prev, y=ETASf_oto_box))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=box_preds), size=1.3, colour="steelblue")+
  ylab("s'")+
  xlab("s")+
  ggtitle("Manual Box-Cox ETASf")+
  theme_classic()



#Extract the residuals from this
ETASf_res_box_size <- residuals(ETASf_box_size_gam)
ETASf_res_box_size <- (ETASf_res_box_size)^2
ETASf <- ETASf %>% mutate(res_box_size=ETASf_res_box_size)

#GAM of the residuals
ETASf_res_box_gam <- gam(ETASf$res_box_size~s(ETASf$prev, k=4))
gam.check(ETASf_res_box_gam)
summary(ETASf_res_box_gam)

ETASf <- ETASf %>% mutate(res_box_preds=predict(ETASf_res_box_gam))
ggplot(ETASf,aes(x=prev, y=res_box_size))+
  geom_point(size=1)+
  geom_line(aes(x=prev, y=res_box_preds), size=1.3, colour="steelblue")+
  ylab("r^2")+
  xlab("s")+
  ggtitle("Manual Box-Cox ETASf")+
  theme_classic()

#Doesn't look any different from the original pre-transformation stuff