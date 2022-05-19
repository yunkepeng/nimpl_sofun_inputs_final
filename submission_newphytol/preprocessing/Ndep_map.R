#output N deposition for 2000-2009.

#now, inputting all predictors
library(tidyverse) 
library(ncmeta)
library(viridis)
library(ggthemes)
library(LSD)
library(yardstick)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gplots)
library(tidyselect)
library(extrafont)
library(caret)
library(recipes)
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")
#library(rbeni)
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)
library(lme4)
library(lmerTest)
library("PerformanceAnalytics")
library(MuMIn)
library(tidyverse)
library(ggplot2)
library(lme4)
library(visreg)
library(ggpubr)
library(car)
library("ggplotify")

vcmax25_df <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vcmax25.nc"),
  varnam = "vcmax25"))

Tg <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/Tg.nc"),
  varnam = "Tg"))

PPFD <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/PPFD.nc"),
  varnam = "PPFD"))

vpd <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vpd.nc"),
  varnam = "vpd"))

fAPAR <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/fAPAR.nc"),
  varnam = "fAPAR"))

age <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/age.nc"),
  varnam = "age"))

CNrt <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/CNrt.nc"),
  varnam = "CNrt"))

LMA <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/LMA.nc"),
  varnam = "LMA"))

###input land cover
ncin <- nc_open("~/data/landcover/modis_landcover_halfdeg_2010_FILLED.nc")
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon) 
lat<-ncvar_get(ncin,"lat")
nlat<-dim(lat)
pftcover <-ncvar_get(ncin,"pftcover")
nc_close(ncin)
pftcover_long <- as.vector(pftcover)
pftcover <- as.data.frame(matrix(pftcover_long, nrow = nlon * nlat, ncol = 10))
#see get_fpc_grid function: https://github.com/stineb/sofun/blob/db7a9e8e486f576fd7b9f1f74edb1df7a8d2c4f7/src/forcing_global_wmodel.mod.f90 
#it clarified that: 1-6 is forest, 8 is grassland
forest_percent <- rowSums(pftcover[,1:6],na.rm=TRUE)/rowSums(pftcover[,c(1:6,8)],na.rm = TRUE)
grass_percent <- pftcover[,8]/rowSums(pftcover[,c(1:6,8)],na.rm = TRUE)
summary(grass_percent + forest_percent) # check - their sum = 1, perfect!

###3. calculate weighted-sum
#firstly - filter na points - so that all output map has same numbers of NA.
all_predictors <- as.data.frame(cbind(Tg$myvar,PPFD$myvar,vpd$myvar,
                                      fAPAR$myvar,age$myvar,
                                      CNrt$myvar,LMA$myvar,vcmax25_df$myvar))
all_predictors$available_grid = rowMeans(all_predictors)
all_predictors$lon <- vcmax25_df$lon
all_predictors$lat <- vcmax25_df$lat
grids <- subset(all_predictors,is.na(available_grid)==FALSE)[,c("lon","lat")]
dim(grids)

devtools::load_all("/Users/yunpeng/yunkepeng/compuetational_ingestr/ingestr/")

grids$nhx <- NA
grids$noy <- NA

for (i in 1:nrow(grids)) {
  tryCatch({
    print(i)
    df_ndep <- ingest_bysite(
      sitename  = paste("a",i,sep=""),
      source    = "ndep",
      lon       = grids$lon[i],
      lat       = grids$lat[i],
      year_start= 2000,
      year_end  = 2009,
      timescale = "y",
      dir       = "~/data/ndep_lamarque/",
      verbose   = FALSE
    )
    grids$noy[i] <- mean(df_ndep$noy,na.rm=TRUE)
    grids$nhx[i] <- mean(df_ndep$nhx,na.rm=TRUE)
  }, error=function(e){})} 

grids$ndep <- grids$nhx+grids$noy

#looks reasonable
plot_map3(grids[,c("lon","lat","ndep")],
          varnam = "ndep",latmin = -65, latmax = 85)


ndep_final <- merge(all_predictors[,c("lon","lat")],grids[,c("lon","lat","ndep")],by=c("lon","lat"),all.x=TRUE)
dim(ndep_final)

library(ncdf4)
ncin <- nc_open("~/data/watch_wfdei/WFDEI-elevation.nc")
lon <- ncvar_get(ncin,"lon")
lat<-ncvar_get(ncin,"lat")

ndep_final_nc <- list(df_to_grid(ndep_final,varnam = "ndep", lonnam = "lon", latnam = "lat"))
names(ndep_final_nc) <- "ndep"
varams = "ndep"
test <- list(lon,lat,ndep_final_nc,varams)
names(test) <- c("lon","lat","vars","varams")
write_nc2(test,varnams = "ndep",long_name = "n_deposition_2000_2009",units = "gN/m2/yr",
          path = "~/data/nimpl_sofun_inputs/map/Final_ncfile/ndep.nc")

d1 <-as.data.frame(nc_to_df(read_nc_onefile("~/data/nimpl_sofun_inputs/map/Final_ncfile/ndep.nc"), varnam = "ndep"))
summary(d1)
plot_map3(na.omit(d1[,c("lon","lat","myvar")]),
          varnam = "myvar",latmin = -65, latmax = 85)
