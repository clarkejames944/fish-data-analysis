source("analysis/setup.R")

library(ggrepel)

###########################################################

# show that the temperatures differ considerably with zone

ylab <- "Seafloor temperature (°C)"
 ggplot(fishdat_cut, aes(x= zone, y= bottomtemp1, fill=zone))+
  geom_boxplot()+
  xlab("Zone") + ylab(ylab)+
  theme_linedraw()
 
  theme_update(legend.position="none")

 

##########################################################    

# plot the body size against otolith size data

fishdat_max <- fishdat_cut %>% group_by(FishID) %>% 
  top_n(1)


ggplot(fishdat_max, aes(x=floorlength, y=z0))+
  geom_point()+
  xlab("Body length at capture (cm)") + ylab("Otolith size at capture(cm)")+
  annotate("text", x=60, y=0.5, label= expression(paste(r^2, "= 0.76, p<0.001")))+
  theme_linedraw()

# obtaining an rsqruared value for the relationship between these measures


cor.test(fishdat_max$floorlength, fishdat_max$z0, method = "pearson")

cor(fishdat_max$floorlength, fishdat_max$z0, method = "pearson")

modlength <- lm(z0~floorlength, data=fishdat_max)

summary.lm(modlength)


#######################################################

### Data summary plots


### bottom temperature against year by zone
fishdat_sum <- fishdat_cut %>% group_by(zone, Year) %>% 
  summarise(bottomtemp1=mean(bottomtemp1))

ylab <- "Seafloor temperature (°C)"
fishdat_sum %>% 
  ggplot(aes(x= Year, y= bottomtemp1, color=zone))+
  xlab("Year") + ylab(ylab)+
  labs(color="Zone")+
  geom_line()



#### Fish ID against year by zone

fishy <- fishdat_cut %>% group_by(zone, FishID) %>% mutate(count=n())
fishdat_arr <- fishdat_cut %>% arrange(zone)
fishdat_arr$Fish.ID <- as.integer(fishdat_arr$FishID)
fishdat_arr$Fish.ID <- as.factor(fishdat_arr$Fish.ID)
fishdat_arr %>% 
  ggplot(aes(x= Year, y=Fish.ID, color=zone))+
  geom_line()
