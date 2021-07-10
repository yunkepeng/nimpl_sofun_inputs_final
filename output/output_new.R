rm(list=ls())
#for figure 1 - just summarizing each part
#load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/output/output.Rdata")

nre_constant_grass <- 0.69

library(tidyverse)  # depends
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
library(rbeni)
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)

firstyr_data <- 1982 # In data file, which is the first year
endyr_data <- 2011 # In data file, which is the last year
location <- "/Users/yunpeng/data/output/latest_forest/"
alloutput_list <- list.files(location,full.names = T)

#input elevation nc file, which will be cbind with global df directly
elev_nc <- read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc")
#elev_nc <- read_nc_onefile("D:/PhD/nimpl_sofun_inputs/Data/Elevation/WFDEI-elevation.nc")
elev <- as.data.frame(nc_to_df(elev_nc, varnam = "elevation"))

#2. Create a function to specify path, loop many years nc file and output a dataframe (lon, lat, var).
inputnc <- function(name,start_year,end_year){
  #-----------------------------------------------------------------------
  # Input: 
  # name: gpp, npp, anpp, vcmax25, leafcn, nuptake...
  # start_year: e.g. 1981
  # end_year: e.g. 2016
  # location: e.g "D:/PhD/nimpl_sofun_inputs/Data/output/" or in Euler: "~/yunkebranch_units/outputnc/"
  #-----------------------------------------------------------------------
  output_allyears <- data.frame(matrix(NA))
  # first, include all years annual data into a daframe
  for (i in firstyr_data:endyr_data){
    if (name == "npp"){
      nc <- read_nc_onefile(alloutput_list[grepl("a.npp.nc", list.files(location,full.names = T))][i-firstyr_data+1]) #we only rely this to filter npp.nc file...
    } else {
      nc <- read_nc_onefile(alloutput_list[grepl(name, list.files(location,full.names = T))][i-firstyr_data+1]) #Input nc
    }
    output_year <- nc_to_df(nc, varnam = name)[,3] #Yearly output
    output_allyears[1:259200,i-firstyr_data+1] <- output_year #here first column represents first year of data file 's output
  }
  names(output_allyears) <- paste(name,firstyr_data:endyr_data,sep="")
  #this variable above (output_allyears), could be end of the function, which is variable at multiple years. But for our purporses, we need mean of select years
  #then, only calculate means of selected years
  output_selected_yrs <- rowMeans(output_allyears[,(start_year-firstyr_data+1):(end_year-firstyr_data+1)],na.rm = TRUE) # only calculated means based on selected start and end year (see function)
  coord <- nc_to_df(nc, varnam = name)[,1:2] # obtain lon and lat
  final_output <- cbind(coord,elev[,3],output_selected_yrs) # combine lon, lat,z with rowmeans variable
  names(final_output) <- c("lon","lat","z",name)
  return(final_output)
  #-----------------------------------------------------------------------
  # Output: output_final: the output data (259200 * 3) including lon, lat and value
  #-----------------------------------------------------------------------
}

#3. select data over 30 years, each df includes lon, lat, z, var
vcmax25_df <- inputnc("vcmax25",1982,2011)

gpp_df <- inputnc("gpp",1982,2011)

#now, inputting all predictors
Tg <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/Tg.nc"),
  varnam = "Tg"))

PPFD <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/PPFD.nc"),
  varnam = "PPFD"))

vpd <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vpd.nc"),
  varnam = "vpd"))

alpha <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/alpha.nc"),
  varnam = "alpha"))

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

#input all regressions 

#firstly, load all forest models
load("/Users/yunpeng/data/NPP_final/statistical_model/mod_tnpp.RData")
summary(mod_tnpp)
load("/Users/yunpeng/data/NPP_final/statistical_model/mod_anpp.RData")
summary(mod_anpp)
load("/Users/yunpeng/data/NPP_final/statistical_model/mod_lnpp.RData")
summary(mod_lnpp)
load("/Users/yunpeng/data/NPP_final/statistical_model/nmass.RData")
summary(n1)
load("/Users/yunpeng/data/NPP_final/statistical_model/nre_model_forest.RData")
summary(nre_model)

#now, do the same for grassland
load(file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/tnpp_grass.RData")
summary(tnpp_grass)
load(file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/anpp_grass.RData")
summary(anpp_grass)

summary(mod_tnpp)$coef
#run global simulations
npp_f <- gpp_df$gpp * (1/(1 + exp(-(summary(mod_tnpp)$coef[1,1]+
                                      summary(mod_tnpp)$coef[2,1]* log(CNrt$myvar)+
                                      summary(mod_tnpp)$coef[3,1] * log(age$myvar) + 
                                      summary(mod_tnpp)$coef[4,1]* fAPAR$myvar))))

anpp_f <- gpp_df$gpp * (1/(1 + exp(-(summary(mod_anpp)$coef[1,1]+
                                       summary(mod_anpp)$coef[2,1] * log(CNrt$myvar)+ 
                                       summary(mod_anpp)$coef[3,1] * log(age$myvar) + 
                                       summary(mod_anpp)$coef[4,1] * fAPAR$myvar))))

bnpp_f <- npp_f-anpp_f

lnpp_f <- anpp_f * (1/(1 + exp(-(summary(mod_lnpp)$coef[1,1]+
                                   summary(mod_lnpp)$coef[2,1]* log(PPFD$myvar) +
                                   summary(mod_lnpp)$coef[3,1] * (Tg$myvar) +
                                   summary(mod_lnpp)$coef[4,1] * log(vpd$myvar)))))

wnpp_f <- anpp_f - lnpp_f

leafcn_f <- (summary(n1)$coef[1,1]/0.46) + 
  (summary(n1)$coef[2,1]/0.46) *vcmax25_df$vcmax25/LMA$myvar

nre_f <- (1/(1+exp(-(summary(nre_model)$coef[1,1]+
                       summary(nre_model)$coef[2,1] *Tg$myvar + 
                       summary(nre_model)$coef[3,1] * log(vpd$myvar)))))

lnf_f <- (1-nre_f)* leafcn_f * lnpp_f

wnf_f <- wnpp_f/100

bnf_f <- bnpp_f/94

nuptake_f <- lnf_f + wnf_f + bnf_f

#grass
npp_g <- gpp_df$gpp * summary(tnpp_grass)$coef[1,1]
anpp_g <- gpp_df$gpp * summary(anpp_grass)$coef[1,1]
bnpp_g <- npp_g-anpp_g
lnf_g <- anpp_g *(1/18)*(1-nre_constant_grass)
bnf_g <- bnpp_g *(1/41)
nuptake_g <- lnf_g + bnf_g



###now, input land cover
ncin <- nc_open("/Users/yunpeng/data/landcover/modis_landcover_halfdeg_2010_FILLED.nc")
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

a1 <- as.data.frame(cbind(gpp_df,forest_percent,grass_percent))

#now, calculate weighted-sum
#firstly - filter na points
all_predictors <- as.data.frame(cbind(Tg$myvar,PPFD$myvar,vpd$myvar,
                                      fAPAR$myvar,age$myvar,
                                      CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25))
all_predictors$available_grid = rowMeans(all_predictors)
#just to find all na columns
all_predictors$available_grid[is.na(all_predictors$available_grid)==FALSE] <- 1
summary(all_predictors$available_grid)
available_grid2 <- all_predictors$available_grid

#represent grids when stand-age is especially in NA, but others are fine
names(all_predictors) <- c("Tg","PPFD","vpd","fAPAR","age","CNrt","LMA","vcmax25","available_grid")
all_predictors$lon <- gpp_df$lon
all_predictors$lat <- gpp_df$lat
summary(all_predictors)


#final calculation
npp_pft <- available_grid2* (npp_f*forest_percent +npp_g*grass_percent)
npp_forest <- available_grid2* (npp_f*forest_percent)
npp_grass <- available_grid2* (npp_g*grass_percent)

anpp_pft <- available_grid2*(anpp_f*forest_percent +anpp_g*grass_percent)
anpp_forest <- available_grid2* (anpp_f*forest_percent)
aanpp_grass <- available_grid2* (anpp_g*grass_percent)

lnpp_forest <- available_grid2*lnpp_f*forest_percent

wnpp_forest <- available_grid2*wnpp_f*forest_percent

bnpp_pft <- available_grid2*(bnpp_f*forest_percent +bnpp_g*grass_percent)
bnpp_forest <- available_grid2* (bnpp_f*forest_percent)
bnpp_grass <- available_grid2* (bnpp_g*grass_percent)

leafcn_forest <- available_grid2*leafcn_f*forest_percent

nre_pft <- available_grid2*nre_f * (forest_percent+grass_percent)

lnf_pft <- available_grid2*(lnf_f*forest_percent +lnf_g*grass_percent)
lnf_forest <- available_grid2* (lnf_f*forest_percent)
lnf_grass <- available_grid2* (lnf_g*grass_percent)

wnf_forest <- available_grid2*wnf_f*forest_percent

bnf_pft <- available_grid2*(bnf_f*forest_percent +bnf_g*grass_percent)
bnf_forest <- available_grid2* (bnf_f*forest_percent)
bnf_grass <- available_grid2* (bnf_g*grass_percent)

nuptake_pft <- available_grid2*(nuptake_f*forest_percent +nuptake_g*grass_percent)
nuptake_pft_final <- nuptake_pft
nuptake_forest <- available_grid2* (nuptake_f*forest_percent)
nuptake_grass <- available_grid2* (nuptake_g*grass_percent)


all_maps <- as.data.frame(cbind(gpp_df,npp_pft,npp_forest,npp_grass,
                                anpp_pft,anpp_forest,aanpp_grass,
                                bnpp_pft,bnpp_forest,bnpp_grass,
                                lnpp_forest,wnpp_forest,leafcn_forest,nre_pft,wnf_forest,
                                lnf_pft,lnf_forest,lnf_grass,
                                bnf_pft,bnf_forest,bnf_grass,
                                nuptake_pft,nuptake_forest,nuptake_grass))
#####area_m2 to show each grid's area in m2
calc_area <- function( lat, dx=1, dy=1 ){
  r_earth <- 6370499.317638  # to be consistent with how Ferret calculates areas of spheres (https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2016/msg00155.html)
  area <- 4 * r_earth^2 * 0.5 * dx * pi/180 * cos( abs(lat) * pi/180 ) * sin( 0.5 * dy * pi/180 )
  return(area)
}
lonlat <- gpp_df[,c("lon","lat")]
area_m2 <- calc_area(lonlat$lat,0.5,0.5)
#fland - to show each grid's land cover percentage
nc <- read_nc_onefile("/Users/yunpeng/data/fland/global.fland.nc") #Input nc
output_fland <- nc_to_df(nc, varnam = "fland")
fland <- output_fland$myvar
#include conversion factor
conversion <- area_m2 * fland /1e+15

#now, defined four things
PAR <- PPFD$myvar*fAPAR$myvar
LUE <- gpp_df$gpp/PAR
CUE <- npp_pft/gpp_df$gpp
NUE_1 <- nuptake_pft/npp_pft
Nup <- PAR*LUE*CUE*NUE_1


plot(NUE_1~CUE)
plot(LUE~NUE_1)

n_PAR <- log(Nup/(mean(PAR,na.rm=TRUE)*LUE*CUE*NUE_1))
n_LUE <- log(Nup/(PAR*mean(LUE,na.rm=TRUE)*CUE*NUE_1))
n_CUE <- log(Nup/(PAR*LUE*mean(CUE,na.rm=TRUE)*NUE_1))
n_NUE_1 <- log(Nup/(PAR*LUE*CUE*mean(NUE_1,na.rm=TRUE)))
n_all <- as.data.frame(cbind(n_PAR,n_LUE,n_CUE,n_NUE_1))
total_sum <- sum(abs((n_all)*conversion),na.rm=TRUE)
#now, research on 4 componenets
sum(abs((n_PAR)*conversion),na.rm=TRUE)/total_sum
sum(abs((n_LUE)*conversion),na.rm=TRUE)/total_sum
sum(abs((n_CUE)*conversion),na.rm=TRUE)/total_sum
sum(abs((n_NUE_1)*conversion),na.rm=TRUE)/total_sum


all_maps2 <- as.data.frame(cbind(gpp_df,PAR,LUE,CUE,NUE_1,Nup))
names(all_maps2)
plot_map3(all_maps2[,c("lon","lat","PAR")],
          varnam = "PAR",
          latmin = -65, latmax = 85,font_size = 15)
ggsave(paste("/Users/yunpeng/data/output/output_map/PAR.jpg",sep=""))
plot_map3(all_maps2[,c("lon","lat","LUE")],
          varnam = "LUE",
          latmin = -65, latmax = 85,font_size = 15)
ggsave(paste("/Users/yunpeng/data/output/output_map/LUE.jpg",sep=""))
plot_map3(all_maps2[,c("lon","lat","CUE")],
          varnam = "CUE",
          latmin = -65, latmax = 85,font_size = 15)
ggsave(paste("/Users/yunpeng/data/output/output_map/CUE.jpg",sep=""))
plot_map3(all_maps2[,c("lon","lat","NUE_1")],
          varnam = "NUE_1",
          latmin = -65, latmax = 85,font_size = 15)
ggsave(paste("/Users/yunpeng/data/output/output_map/NUE_1.jpg",sep=""))
plot_map3(all_maps2[,c("lon","lat","Nup")],
          varnam = "Nup",
          latmin = -65, latmax = 85,font_size = 15)
ggsave(paste("/Users/yunpeng/data/output/output_map/Nup.jpg",sep=""))

My_Theme = theme(
  axis.title.x = element_text(size = 40),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 40),
  axis.text.y = element_text(size = 20))

ggplot(data=all_maps2, aes(x=LUE, y=CUE)) +
  geom_point()+geom_smooth(method = "lm", se = TRUE)+
  xlab("LUE")+ylab("CUE")+theme_classic()+My_Theme
summary(lm(CUE~LUE,all_maps2))

ggplot(data=all_maps2, aes(x=LUE, y=NUE_1)) +
  geom_point()+geom_smooth(method = "lm", se = TRUE)+
  xlab("LUE")+ylab("1 / NUE")+theme_classic()+My_Theme
summary(lm(NUE_1~LUE,all_maps2))

ggplot(data=all_maps2, aes(x=CUE, y=NUE_1)) +
  geom_point()+geom_smooth(method = "lm", se = TRUE)+
  xlab("CUE")+ylab("1 / NUE")+theme_classic()+My_Theme
summary(lm(NUE_1~CUE,all_maps2))
