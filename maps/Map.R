#######################################################################
## Load and attach libraries ----

library(dplyr)
library(ggplot2)
library(rstan)
library(ggmcmc)
library(bayesplot)
library(scales)
library(coda)
library(stringr)
library(mgcv)
library(loo)
library(car)
library(Matrix)
library(viridis)
library(reshape2)
library(tidyverse)
library(maps)
library(mapdata)
library(sp)
library(oz)
library(oceanmap)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggmap)
library(ggsn)
library(ggspatial)
library(anchors)
########################################################################
outline <- maps::map("worldHires", regions=c("Australia", "Tasmania"), exact=TRUE, plot=FALSE) # returns a list of x/y coords
xrange <- range(outline$x, na.rm=TRUE) # get bounding box
yrange <- range(outline$y, na.rm=TRUE)
xbox <- xrange + c(-2, 2)
ybox <- yrange + c(-2, 2)

subset <- !is.na(outline$x)

polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="light blue", rule="evenodd")

########################################################################
oz(states = TRUE, coast = TRUE, xlim = NULL,
   ylim = NULL, add = FALSE, ar = 1, eps = 0.25,
   sections = NULL, visible = NULL)


aus <- oz(states = FALSE, coast = TRUE, xlim = NULL,
   ylim = NULL, add = FALSE, ar = 1, eps = 0.25,
   sections = c(4,5,6,7,8,10, 12,13,14,15,16), visible = NULL)

oz(sections=c(4,5,6,7,8,10, 12,13,14,15,16))

########################################################################
aus <- maps::map("worldHires", "Australia", fill=TRUE, xlim=c(110,160),
         ylim=c(-45,-5), mar=c(0,0,0,0))
aus <- fortify(aus)
ins_map <- ggplot() + 
  geom_polygon(data= aus, aes(x=long, y=lat, group=group))+
  coord_equal()+
  theme_classic()+labs(x=NULL, y=NULL)+
  geom_rect(data=pol, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0, size=1, linetype=1)


south_aus <- maps::map("worldHires", "Australia", fill=TRUE, xlim=c(130,160),
               ylim=c(-45,-20), mar=c(0,0,0,0))

south_aus <- fortify(south_aus)

points(x=map_points$centlong, y=map_points$centlat, col="red", pch=20)
text(x=map_points$centlong, y=map_points$centlat, map_points$zone, pos=4)

map <- ggplot(south_aus, aes(long, lat)) + 
  geom_polygon(aes(group=group), fill="khaki")+
  coord_equal()+
  scale_y_continuous(limits= c(-45,-32.5), expand = c( 0 , 0 ))+
  scale_x_continuous(limits=c(130,155), expand = c( 0 , 0 ))+
  geom_point(data=map_points, color="steelblue")+
  geom_text(data=map_points, label=map_points$zone, vjust=1.5)+
  ggsn::scalebar(south_aus, dist=400, dist_unit = "km", transform = TRUE, height=0.0015, box.fill = c("black", "white"), model = "International")+
  theme_classic()

  ggsn::north2(map, x = 0.1, y=0.3)


  pol<-data.frame(xmin=5,xmax=10 ,ymin=-1 ,ymax=4)
  
  grid.newpage()
  v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
  v2<-viewport(width = 0.3, height = 0.3, x = 0.86, y = 0.28) #plot area for the inset map
  print(map,vp=v1) 
  print(ins_map,vp=v2)


#######################################################################

#Create table for centre of latitiude and longitude for each zone

map_points <- fishdat_cut %>% dplyr::select(zone, centlat, centlong)
map_points <- map_points %>% group_by(zone) %>% 
  summarise(lat=min(centlat),
            long=min(centlong))
map_points[2,3] <- 148.5
map_points[4,3] <- 145


########################################################################
aus <- ne_countries(country = "australia", scale = 'large', returnclass = 'sf')

uk_sac <- "SESS/SESS.shp" %>% 
  st_read() 

aus_sac <- "SESS Commonwealth Trawl-2/SESS Commonwealth Trawl.shp" %>% 
  st_read()

uk_sac <- "csq_tw_20112014/csq_tw_20112014.shp" %>% 
  st_read() 

uk_sac
aus_sac

ggplot(uk_sac) + geom_sf()

st_crs(uk_sac)

st_crs(aus)

ggplot(aus) + geom_sf(fill='grey10')+
  geom_sf(data=uk_sac, fill="orange", alpha=0.5)

ggplot(aus) + geom_sf(fill='grey10')+
  geom_sf(data=aus_sac, fill="orange", alpha=0.5)

######################################################################

