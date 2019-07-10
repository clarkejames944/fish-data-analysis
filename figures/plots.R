source("analysis/setup.R")

###########################################################

# show that the temperatures differ considerably with zone

ggplot(fishdat_cut, aes(x= bottomtemp1, y= z0, color=zone))+
         geom_point()
ylab <- "Seafloor temperature (°C)"
ggplot(fishdat_cut, aes(x= zone, y= bottomtemp1))+
  geom_boxplot()+
  xlab("Zone") + ylab(ylab)+
  theme_classic()

##########################################################    

# plot the body size against otolith size data

fishdat_max <- fishdat_cut %>% group_by(FishID) %>% 
  top_n(1)

fishdat_max <- fishdat_max %>% filter(zone=="EBS", sex=="F")

ggplot(fishdat_max, aes(x=floorlength, y=z0))+
  geom_point()+
  xlab("Body length (cm)") + ylab("Otolith size (cm)")+
  theme_classic()

# obtaining an rsqruared value for the relationship between these measures

modlength <- lm(z0~floorlength, data=fishdat_max)

summary.lm(modlength)
