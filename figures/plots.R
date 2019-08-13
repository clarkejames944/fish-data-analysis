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
fishdat_cut$Year <- as.numeric(fishdat_cut$Year)

fishdat_arr <- fishdat_cut %>% group_by(FishID) %>% arrange(min(Year))

Fish.ID <- 1:length(unique(fishdat_arr$FishID))

fishdat_arr$Fish.ID <- Fish.ID[fishdat_arr$FishID]

test <- NULL
x <- split(fishdat_cut, fishdat_cut$zone)
for(i in seq(x)){
  if(i == 1){
    s <- x[[i]]
    s$ID <- as.numeric(as.factor(as.character(s$FishID)))
    s <- s[order(s$Year),]
    test <- c(test, unique(s$ID))
    
    x[[i]] <- s
  } else {
    s <- x[[i]]
    
    s$ID <- (as.numeric(as.factor(as.character(s$FishID)))) + max(x[[(i-1)]]$ID)
    s <- s[order(s$Year),]
    test <- c(test, unique(s$ID))
    
    x[[i]] <- s
  }
}
x <- do.call("rbind", x)
lapply(split(x, x$zone), function(x) c(head(x$ID), tail(x$ID)))


x$ID <- factor(x$ID, levels = test)
levels(x$ID)
test

x$ID <- as.factor(x$ID)
# x$ID <- droplevels(x$ID)
# x$ID <- factor(x$ID, levels = unique(names(sort(sapply(split(x, x$ID), function(x) min(x$Year))))))

ggplot(data = x, aes(x= Year, y=ID, color=zone))+
  geom_line() +
  geom_point(size = 0.1)+
  # facet_wrap("zone", scales = "free_y")+
  scale_y_discrete(breaks = seq(0, 5000 , by = 200))+
  geom_text(x=1975, y=1500, label="EBS")+
  geom_text(x=1975, y=3700, label="ETAS")+
  geom_text(x=1975, y=3300, label="WTAS")+
  geom_text(x=1975, y=4000, label="NSW")+
  theme(axis.text.y=element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())


