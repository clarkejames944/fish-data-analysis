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
aus<-maps::map("worldHires", "Australia", fill=TRUE, xlim=c(110,160),
         ylim=c(-45,-5), mar=c(0,0,0,0))

ggplot(fortify(aus), aes(y=lat, x=long, group=group)) + geom_polygon()


south_aus<-maps::map("worldHires", "Australia", fill=TRUE, xlim=c(130,160),
               ylim=c(-45,-20), mar=c(0,0,0,0))

ggplot(fortify(south_aus), aes(y=lat, x=long,)) + geom_polygon()

aus.sp <- maps::map2SpatialPolygons(aus)
par(mar=c(0,0,0,0))
plot(aus.sp, asp=1)
ggplot(fortify(aus.sp), aes(y=lat, x=long, group=group)) + geom_polygon()

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
