rm(list=ls())

load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/output/output.Rdata")

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
#forest map
plot_map3(a1[,c("lon","lat","forest_percent")],
          varnam = "forest_percent",
          latmin = -65, latmax = 85)
#grass map
plot_map3(a1[,c("lon","lat","grass_percent")],
          varnam = "grass_percent",
          latmin = -65, latmax = 85)

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
npp_pft_final <- npp_pft
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

#just print values
for (i in 4:ncol(all_maps)){
  varname <- names(all_maps)[i]
  total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
  print(varname)
  print(total_value)
}


#plot Nuptake ~ gpp
My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 20))

simulations <- all_maps
simulations$nue <- simulations$nuptake_pft/simulations$npp_pft
#simulations$nue_gpp <- simulations$nuptake_pft/simulations$gpp
simulations_predictors <- as.data.frame(cbind(all_predictors[,-c(9,10,11)],simulations))

#####nue
plot_map3(simulations[,c("lon","lat","nuptake_pft")],
          varnam = "nuptake_pft",
          latmin = -65, latmax = 85)
plot_map3(simulations[,c("lon","lat","nuptake_forest")],
          varnam = "nuptake_pft",
          latmin = -65, latmax = 85)

plot_map3(simulations[,c("lon","lat","nue")],
          varnam = "nue",
          latmin = -65, latmax = 85)

#####nue -forest, grassland
simulations$nue_forest <- simulations$nuptake_forest/simulations$npp_forest
simulations$nue_grass <- simulations$nuptake_grass/simulations$npp_grass

plot_map3(simulations[,c("lon","lat","nue_forest")],
          varnam = "nue_forest",
          latmin = -65, latmax = 85)

ggplot(data=simulations, aes(simulations$nue)) + 
  geom_histogram()+
  geom_vline(aes(xintercept = mean(simulations$nue,na.rm=TRUE),col='all'),size=2)+
  geom_vline(aes(xintercept = mean(simulations$nuptake_forest/simulations$npp_forest,na.rm=TRUE),col='forest'),size=2)+
  geom_vline(aes(xintercept = mean(simulations$nuptake_grass/simulations$npp_grass,na.rm=TRUE),col='grassland'),size=2)+
  theme(text = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))+
  theme(legend.title = element_blank())+
  labs(x = "N demand of NPP (gN/gC)")

hist(1/simulations$nue_gpp)
#hist(1/simulations$nue_gpp_grass)
#hist(1/simulations$nue_gpp_forest)

ggplot(data=simulations, aes(nuptake_pft)) + 
  geom_histogram()+
  geom_vline(aes(xintercept = mean(nuptake_pft,na.rm=TRUE),col='all'),size=2)+
  geom_vline(aes(xintercept = mean(nuptake_f,na.rm=TRUE),col='forest'),size=2)+
  geom_vline(aes(xintercept = mean(nuptake_g,na.rm=TRUE),col='grassland'),size=2)+
  theme(text = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))+
  theme(legend.title = element_blank())+
  labs(x = "Total N uptake (gN/m2/yr)") 


#now, select median values from available grids only - and do each step-by-step
length(Tg$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_Tg <- mean(Tg$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_PPFD <- mean(PPFD$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_vpd <- mean(vpd$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_fAPAR <- mean(fAPAR$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_age <- mean(age$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_CNrt <- mean(CNrt$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_LMA <- mean(LMA$myvar[is.na(all_predictors$available_grid)==FALSE])
mean_vcmax25 <- mean(vcmax25_df$vcmax25[is.na(all_predictors$available_grid)==FALSE])

#create a function to control each factor --> output ratio of N demand of NPP
cal_nuptake <- function(Tg_pred,PPFD_pred,vpd_pred,fAPAR_pred,age_pred,CNrt_pred,LMA_pred,vcmax25_pred){
  npp_f <- gpp_df$gpp * (1/(1 + exp(-(summary(mod_tnpp)$coef[1,1]+
                                        summary(mod_tnpp)$coef[2,1]* log(CNrt_pred)+
                                        summary(mod_tnpp)$coef[3,1] * log(age_pred) + 
                                        summary(mod_tnpp)$coef[4,1]* fAPAR_pred))))
  
  anpp_f <- gpp_df$gpp * (1/(1 + exp(-(summary(mod_anpp)$coef[1,1]+
                                         summary(mod_anpp)$coef[2,1] * log(CNrt_pred)+ 
                                         summary(mod_anpp)$coef[3,1] * log(age_pred) + 
                                         summary(mod_anpp)$coef[4,1] * fAPAR_pred))))
  
  bnpp_f <- npp_f-anpp_f
  
  lnpp_f <- anpp_f * (1/(1 + exp(-(summary(mod_lnpp)$coef[1,1]+
                                     summary(mod_lnpp)$coef[2,1]* log(PPFD_pred) +
                                     summary(mod_lnpp)$coef[3,1] * (Tg_pred) +
                                     summary(mod_lnpp)$coef[4,1] * log(vpd_pred)))))
  
  wnpp_f <- anpp_f - lnpp_f
  
  leafcn_f <- (summary(n1)$coef[1,1]/0.46) + 
    (summary(n1)$coef[2,1]/0.46) *vcmax25_pred/LMA_pred
  
  nre_f <- (1/(1+exp(-(summary(nre_model)$coef[1,1]+
                         summary(nre_model)$coef[2,1] *Tg_pred + 
                         summary(nre_model)$coef[3,1] * log(vpd_pred)))))
  
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
  
  #pft
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
  nuptake_forest <- available_grid2* (nuptake_f*forest_percent)
  nuptake_grass <- available_grid2* (nuptake_g*grass_percent)
  nuptake_pft_ratio <- (nuptake_pft_final/npp_pft_final)/(nuptake_pft/npp_pft)
  return(nuptake_pft_ratio)
}

nuptake_standard <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
summary(nuptake_standard)

nuptake_Tg <- cal_nuptake(rep(mean_Tg,259200),PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_PPFD <- cal_nuptake(Tg$myvar,rep(mean_PPFD,259200),vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_vpd <- cal_nuptake(Tg$myvar,PPFD$myvar,rep(mean_vpd,259200),fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_fAPAR <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,rep(mean_fAPAR,259200),age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_age <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,rep(mean_age,259200),CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_CNrt <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,rep(mean_CNrt,259200),LMA$myvar,vcmax25_df$vcmax25)
nuptake_LMA <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,rep(mean_LMA,259200),vcmax25_df$vcmax25)
nuptake_vcmax25 <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,rep(mean_vcmax25,259200))

nuptake_all <- as.data.frame(cbind(nuptake_Tg,nuptake_PPFD,nuptake_vpd,nuptake_fAPAR,nuptake_age,nuptake_CNrt,nuptake_LMA,nuptake_vcmax25))
summary(nuptake_all)

#output most max factor, and its value
for (i in 1:nrow(nuptake_all)){
  if (is.na(nuptake_all[i,1])==TRUE){
    nuptake_all$most_factor[i] <- NA
    nuptake_all$most_factor_value[i] <- NA
  } else {
    nuptake_all$most_factor[i] <- names((which.max(abs(log(nuptake_all[i,1:8])))))
    nuptake_all$most_factor_value[i] <- max(abs(log(nuptake_all[i,1:8])))
  }
}

apparent_point <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all[,c("most_factor","most_factor_value")]))
apparent_point_available <- subset(apparent_point,is.na(most_factor_value)==FALSE) 
dim(apparent_point_available)
apparent_point_available$most_factor_value <- as.numeric(apparent_point_available$most_factor_value)

apparent_point_available %>% group_by(most_factor) %>% summarise(number=n())

area_final <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all[,c("most_factor","most_factor_value")]))
area_final$conversion <- conversion

dim(subset(area_final,is.na(most_factor_value)==FALSE)) # area only available in those grids
#total area
sum(subset(area_final,is.na(most_factor_value)==FALSE)$conversion, na.rm = TRUE)

area_final$ratio <- area_final$conversion/sum(subset(area_final,is.na(most_factor_value)==FALSE)$conversion, na.rm = TRUE)

aa <- area_final %>% group_by(most_factor) %>% summarise(sum = sum(ratio)*100, n = n())
sum(aa$sum,na.rm=TRUE)
aa

# count the needed levels of a factor
gg <- plot_map3(all_maps[,c("lon","lat","nuptake_pft")],
                varnam = "nuptake_pft",plot_title = paste("N demand of NPP constrained by most important factor"),
                latmin = -65, latmax = 85,combine=FALSE)

apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_age"] <- "age"
apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_CNrt"] <- "Soil C/N"
apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_fAPAR"] <- "fAPAR"
apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_LMA"] <- "LMA"
apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_PPFD"] <- "PPFD"
apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_Tg"] <- "Tg"
apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_vcmax25"] <- "Vcmax25"
apparent_point_available$most_factor[apparent_point_available$most_factor=="nuptake_vpd"] <- "vpd"

colors <-  c("red","red","red","red","red","red","red","red","cyan","black","yellow","green","orange","red","blue","purple")
gg$ggmap + geom_point(data=apparent_point_available,aes(lon,lat,color=most_factor),size=0.5)+
  scale_color_manual(values = colors)+ theme(
    legend.text = element_text(size = 20))+
    guides(colour = guide_legend(override.aes = list(size = 5)))
ggsave(paste("/Users/yunpeng/data/output/output_onefactor/allfactor.jpg",sep=""))

#now, output each factor's effect
#devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/rbeni/")

nuptake_all_coord <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all))
summary(nuptake_all_coord)
nuptake_all_coord[,3:10] <- log(nuptake_all_coord[,3:10])
summary(nuptake_all_coord[,3:10])

total_sum <- sum(abs((nuptake_all_coord[,3:10])*conversion),na.rm=TRUE)
total_sum

#global trend
global_trend <- (nuptake_all_coord[,3:10])*conversion/sum(conversion,na.rm = TRUE)
colSums(global_trend,na.rm=TRUE)
sum(colSums(global_trend,na.rm=TRUE))
exp(sum(colSums(global_trend,na.rm=TRUE)))
#increased 6%

names(nuptake_all_coord) <- c("lon","lat","Tg","PPFD","vpd","fAPAR","age","Soil C:N","LMA","Vcmax25","most_factor","most_factor_value")

for (i in c(3,5,10)){
  varname <- names(nuptake_all_coord)[i]
  relative_value <- round(sum(abs((nuptake_all_coord[,i])*conversion),na.rm=TRUE)/total_sum,2)
  percentage_value <- label_percent()(relative_value)
  plot_map3(nuptake_all_coord[,c("lon","lat",varname)],
            varnam = varname,plot_title = paste(varname, percentage_value, sep=": " ),
            latmin = -65, latmax = 85,
            colorscale = c( "royalblue4", "wheat","tomato3"),
            breaks = c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.025, 0,
                       0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5))
  ggsave(paste("/Users/yunpeng/data/output/output_onefactor/",varname,".jpg",sep=""))
}

for (i in c(4,6,7,8,9)){
  varname <- names(nuptake_all_coord)[i]
  relative_value <- round(sum(abs((nuptake_all_coord[,i])*conversion),na.rm=TRUE)/total_sum,2)
  percentage_value <- label_percent()(relative_value)
  plot_map3(nuptake_all_coord[,c("lon","lat",varname)],
            varnam = varname,plot_title = paste(varname, percentage_value, sep=": " ),
            latmin = -65, latmax = 85,
            colorscale = c( "royalblue4", "wheat","tomato3"),
            breaks = c(-0.1,-0.05,-0.025, 0,
                       0.025,0.05,0.1))
  ggsave(paste("/Users/yunpeng/data/output/output_onefactor/",varname,".jpg",sep=""))
}

save.image(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/output/output.Rdata")


######
#quantify each component
#create a function to control each factor --> output ratio of nuptake_factor / nuptake_standard

model1_npp_gpp <- available_grid2*(1/(1 + exp(-(summary(mod_tnpp)$coef[1,1]+
                                                  summary(mod_tnpp)$coef[2,1]* log(CNrt$myvar)+
                                                  summary(mod_tnpp)$coef[3,1] * log(age$myvar) + 
                                                  summary(mod_tnpp)$coef[4,1]* fAPAR$myvar))))

model2_anpp_gpp <- available_grid2*(1/(1 + exp(-(summary(mod_anpp)$coef[1,1]+
                                                   summary(mod_anpp)$coef[2,1] * log(CNrt$myvar)+ 
                                                   summary(mod_anpp)$coef[3,1] * log(age$myvar) + 
                                                   summary(mod_anpp)$coef[4,1] * fAPAR$myvar))))

model3_lnpp_anpp <- available_grid2*(1/(1 + exp(-(summary(mod_lnpp)$coef[1,1]+
                                                    summary(mod_lnpp)$coef[2,1]* log(PPFD$myvar) +
                                                    summary(mod_lnpp)$coef[3,1] * (Tg$myvar) +
                                                    summary(mod_lnpp)$coef[4,1] * log(vpd$myvar)))))

model4_leafnc <- available_grid2*(summary(n1)$coef[1,1]/0.46) + 
  (summary(n1)$coef[2,1]/0.46) *vcmax25_df$vcmax25/LMA$myvar

model5_nre <- available_grid2*(1/(1+exp(-(summary(nre_model)$coef[1,1]+
                                            summary(nre_model)$coef[2,1] *Tg$myvar + 
                                            summary(nre_model)$coef[3,1] * log(vpd$myvar)))))

cal_nuptake_model <- function(gpp,model1_npp_gpp,model2_anpp_gpp,model3_lnpp_anpp,model4_leafnc,model5_nre){
  npp_f <- gpp* model1_npp_gpp
  
  anpp_f <- gpp * model2_anpp_gpp
  bnpp_f <- npp_f-anpp_f
  
  lnpp_f <- anpp_f * model3_lnpp_anpp  
  wnpp_f <- anpp_f - lnpp_f
  
  leafcn_f <- model4_leafnc
  
  nre_f <- model5_nre 
  lnf_f <- (1-nre_f)* leafcn_f * lnpp_f
  
  wnf_f <- wnpp_f/100
  
  bnf_f <- bnpp_f/94
  
  nuptake_f <- lnf_f + wnf_f + bnf_f
  
  #grass
  npp_g <- gpp* summary(tnpp_grass)$coef[1,1]
  anpp_g <- gpp * summary(anpp_grass)$coef[1,1]
  bnpp_g <- npp_g-anpp_g
  lnf_g <- anpp_g *(1/18)*(1-nre_constant_grass)
  bnf_g <- bnpp_g *(1/41)
  nuptake_g <- lnf_g + bnf_g
  
  #pft
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
  nuptake_forest <- available_grid2* (nuptake_f*forest_percent)
  nuptake_grass <- available_grid2* (nuptake_g*grass_percent)
  nuptake_pft_ratio <- nuptake_pft_final/nuptake_pft
  return(nuptake_pft_ratio)
}

nuptake_standard <- cal_nuptake_model(gpp_df$gpp,model1_npp_gpp,model2_anpp_gpp,model3_lnpp_anpp,model4_leafnc,model5_nre)
summary(nuptake_standard)

mean_gpp <- mean(gpp_df$gpp[is.na(all_predictors$available_grid)==FALSE])
mean_npp_gpp <- mean(model1_npp_gpp[is.na(all_predictors$available_grid)==FALSE])
mean_anpp_gpp <- mean(model2_anpp_gpp[is.na(all_predictors$available_grid)==FALSE])
mean_lnpp_anpp <- mean(model3_lnpp_anpp[is.na(all_predictors$available_grid)==FALSE])
mean_leafnc <- mean(model4_leafnc[is.na(all_predictors$available_grid)==FALSE])
mean_nre <- mean(model5_nre[is.na(all_predictors$available_grid)==FALSE])

gpp_model <- cal_nuptake_model(rep(mean_gpp,259200),model1_npp_gpp,model2_anpp_gpp,model3_lnpp_anpp,model4_leafnc,model5_nre)
npp_gpp_model <- cal_nuptake_model(gpp_df$gpp,rep(mean_npp_gpp,259200),model2_anpp_gpp,model3_lnpp_anpp,model4_leafnc,model5_nre)
anpp_gpp_model <- cal_nuptake_model(gpp_df$gpp,model1_npp_gpp,rep(mean_anpp_gpp,259200),model3_lnpp_anpp,model4_leafnc,model5_nre)
lnpp_anpp_model <- cal_nuptake_model(gpp_df$gpp,model1_npp_gpp,model2_anpp_gpp,rep(mean_lnpp_anpp,259200),model4_leafnc,model5_nre)
leafnc_model <- cal_nuptake_model(gpp_df$gpp,model1_npp_gpp,model2_anpp_gpp,model3_lnpp_anpp,rep(mean_leafnc,259200),model5_nre)
nre_model <- cal_nuptake_model(gpp_df$gpp,model1_npp_gpp,model2_anpp_gpp,model3_lnpp_anpp,model4_leafnc,rep(mean_nre,259200))

nuptake_all2 <- as.data.frame(cbind(npp_gpp_model,anpp_gpp_model,
                                    lnpp_anpp_model,leafnc_model,nre_model))
#disregard gpp effect - now it researched on how each model affect Nup/gpp
summary(nuptake_all2)

for (i in 1:nrow(nuptake_all2)){
  if (is.na(nuptake_all2[i,1])==TRUE){
    nuptake_all2$most_factor[i] <- NA
    nuptake_all2$most_factor_value[i] <- NA
  } else { 
    nuptake_all2$most_factor[i] <- names((which.max(abs(log(nuptake_all2[i,1:5])))))
    nuptake_all2$most_factor_value[i] <- max(abs(log(nuptake_all2[i,1:5])))
  }
}

apparent_point2 <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all2[,c("most_factor","most_factor_value")]))
apparent_point_available2 <- subset(apparent_point2,is.na(most_factor_value)==FALSE) 
dim(apparent_point_available2)

apparent_point_available2 %>% group_by(most_factor) %>% dplyr::summarise(number=n())

area_final2 <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all2[,c("most_factor","most_factor_value")]))
area_final2$conversion <- conversion

sum(subset(area_final2,is.na(most_factor_value)==FALSE)$conversion, na.rm = TRUE)

area_final2$ratio <- area_final2$conversion/sum(subset(area_final2,is.na(most_factor_value)==FALSE)$conversion, na.rm = TRUE)

bb <- area_final2 %>% group_by(most_factor) %>% summarise(sum = sum(ratio)*100, n = n())
sum(bb$sum,na.rm=TRUE)
bb

# count the needed levels of a factor
gg <- plot_map3(all_maps[,c("lon","lat","nuptake_pft")],
                varnam = "nuptake_pft",plot_title = paste("N uptake/GPP constrained by most important model"),
                latmin = -65, latmax = 85,combine=FALSE)

colors <-  c("red","red","red","red","red","red","red","red", "blue","yellow","red","black","orange","red","green","purple")
gg$ggmap + geom_point(data=apparent_point_available2,aes(lon,lat,color=most_factor),size=0.5)+
  scale_color_manual(values = colors)+ theme(
    legend.text = element_text(size = 20))+
  guides(colour = guide_legend(override.aes = list(size = 5)))
ggsave(paste("/Users/yunpeng/data/output/output_onefactor/allfactor_model.jpg",sep=""))

nuptake_all_coord2 <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all2))
summary(nuptake_all_coord2)
nuptake_all_coord2[,3:8] <- log(nuptake_all_coord2[,3:8])
summary(nuptake_all_coord2[,3:8])

#now, output each factor
total_sum2 <- sum(abs((nuptake_all_coord2[,3:8])*conversion),na.rm=TRUE)
total_sum2

#global_trend
global_trend2 <- (nuptake_all_coord2[,3:8])*conversion/sum(conversion,na.rm = TRUE)
colSums(global_trend2,na.rm=TRUE)
sum(colSums(global_trend2,na.rm=TRUE))
exp(sum(colSums(global_trend2,na.rm=TRUE)))
#increased 6%


for (i in 3:8){
  varname <- names(nuptake_all_coord2)[i]
  relative_value <- round(sum(abs((nuptake_all_coord2[,i])*conversion),na.rm=TRUE)/total_sum2,2)
  percentage_value <- label_percent()(relative_value)
  plot_map3(nuptake_all_coord2[,c("lon","lat",varname)],
            varnam = varname,plot_title = paste(varname, percentage_value, sep=": " ),
            latmin = -65, latmax = 85,breaks=c(seq(-4, 2, by = 0.5)),
            colorscale = c( "royalblue4","royalblue2", "wheat", "tomato2" ))
  ggsave(paste("/Users/yunpeng/data/output/output_onefactor/",varname,".jpg",sep=""))
}
#save.image(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/output/output.Rdata")


#now, deal with how npp/gpp...nre would affect nuptake/gpp
nuptake_all3 <- as.data.frame(cbind(npp_gpp_model,anpp_gpp_model,
                                    lnpp_anpp_model,leafnc_model,nre_model))

for (i in 1:nrow(nuptake_all3)){
  if (is.na(nuptake_all3[i,1])==TRUE){
    nuptake_all3$most_factor[i] <- NA
    nuptake_all3$most_factor_value[i] <- NA
  } else { 
    nuptake_all3$most_factor[i] <- names((which.max(abs(log(nuptake_all3[i,1:5])))))
    nuptake_all3$most_factor_value[i] <- max(abs(log(nuptake_all3[i,1:5])))
  }
}

apparent_point3 <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all3[,c("most_factor","most_factor_value")]))
apparent_point_available3 <- subset(apparent_point3,is.na(most_factor_value)==FALSE) 
dim(apparent_point_available2)

apparent_point_available3 %>% group_by(most_factor) %>% summarise(number=n())

area_final3 <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all3[,c("most_factor","most_factor_value")]))
area_final3$conversion <- conversion

sum(subset(area_final3,is.na(most_factor_value)==FALSE)$conversion, na.rm = TRUE)

area_final3$ratio <- area_final3$conversion/sum(subset(area_final3,is.na(most_factor_value)==FALSE)$conversion, na.rm = TRUE)

cc <- area_final3 %>% group_by(most_factor) %>% summarise(sum = sum(ratio)*100, n = n())
sum(cc$sum,na.rm=TRUE)
cc

# count the needed levels of a factor
gg <- plot_map3(all_maps[,c("lon","lat","nuptake_pft")],
                varnam = "nuptake_pft",plot_title = paste("N uptake Constrained by most important factor"),
                latmin = -65, latmax = 85,combine=FALSE)

colors <-  c("red","red","red","red","red","red","red","red", "brown","royalblue3","plum2","red","yellow","red","red","red")
gg$ggmap + geom_point(data=apparent_point_available3,aes(lon,lat,color=most_factor),size=0.5)+
  scale_color_manual(values = colors)+ theme(
    legend.text = element_text(size = 20))+
  guides(colour = guide_legend(override.aes = list(size = 5)))
ggsave(paste("/Users/yunpeng/data/output/output_onefactor/allfactor_model_outgpp.jpg",sep=""))

nuptake_all_coord3 <- as.data.frame(cbind(gpp_df[,c("lon","lat")],nuptake_all3))
summary(nuptake_all_coord3)
nuptake_all_coord3[,3:7] <- log(nuptake_all_coord3[,3:7])
summary(nuptake_all_coord3[,3:7])

#now, output each factor
total_sum3 <- sum(abs((nuptake_all_coord3[,3:7])*conversion),na.rm=TRUE)
total_sum3

for (i in 3:7){
  varname <- names(nuptake_all_coord3)[i]
  relative_value <- round(sum(abs((nuptake_all_coord3[,i])*conversion),na.rm=TRUE)/total_sum3,2)
  percentage_value <- label_percent()(relative_value)
  plot_map3(nuptake_all_coord3[,c("lon","lat",varname)],
            varnam = varname,plot_title = paste(varname, percentage_value, sep=": " ),
            latmin = -65, latmax = 85,
            colorscale = c( "royalblue4", "wheat", "tomato2"),
            breaks = c(-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-0.025, 0,
                       0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45))
  ggsave(paste("/Users/yunpeng/data/output/output_onefactor/","more_",varname,".jpg",sep=""))
}

#now, deal with uncertainty
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

## Uncertainty of TNPP/GPP
#using standard error method to calculate uncertainty of whole regression
# for all logit function model (mod_tnpp, mod_anpp, mod_lnpp, nre_model): a = 1/(1+exp(-b)), where a is the ratio (e.g. npp/gpp) and b is the regression (e.g. mod =tnpp)
#we calculate uncertainty of b firstly, based on each regression
# npp/gpp uncertainty
mod_tnpp_uncertainty <- summary(mod_tnpp)$sigma  #random factor: a population parameter - estimated automatically
#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
tnpp_b <- summary(mod_tnpp)$coefficients[1,1] + summary(mod_tnpp)$coefficients[2,1] * log(CNrt$myvar)  +
  summary(mod_tnpp)$coefficients[3,1] * log(age$myvar) + summary(mod_tnpp)$coefficients[4,1] * fAPAR$myvar
#therefore, the uncertainty of npp / gpp (defined as a here) = uncertainty b * (deriative a / deriative b)
mod_tnpp_uncertainty <- mod_tnpp_uncertainty *  exp(-tnpp_b) / ( (1 + exp(-tnpp_b)) ^2)
summary(mod_tnpp_uncertainty)

## Uncertainty of ANPP/GPP
mod_anpp_uncertainty <- summary(mod_anpp)$sigma
anpp_b <- summary(mod_anpp)$coefficients[1,1] + summary(mod_anpp)$coefficients[2,1] * log(CNrt$myvar)  +
  summary(mod_anpp)$coefficients[3,1] * log(age$myvar) + summary(mod_anpp)$coefficients[4,1] * fAPAR$myvar
mod_anpp_uncertainty <- mod_anpp_uncertainty *  exp(-anpp_b) / ( (1 + exp(-anpp_b)) ^2)


## Uncertainty of LNPP/ANPP
mod_lnpp_uncertainty <- summary(mod_lnpp)$sigma
lnpp_b <- summary(mod_lnpp)$coefficients[1,1] + summary(mod_lnpp)$coefficients[2,1] * log(PPFD$myvar)  +
  summary(mod_lnpp)$coefficients[3,1] * Tg$myvar + summary(mod_lnpp)$coefficients[4,1] * log(vpd$myvar)
mod_lnpp_uncertainty <- mod_lnpp_uncertainty *  exp(-lnpp_b) / ( (1 + exp(-lnpp_b)) ^2)

#leaf n/c - ALREADY checked that inputted vcmax25_df + LMA calculated from (1) fortran and (2) R for predicting leaf n/c, could output the same prediction for leaf n/c Please note! using a.vcmax25.nc rather than annualvcmax25.nc. See difference in: https://www.notion.so/computationales/annualvcmax25-vs-a-vcmax25-282bddd63bba4d7b9cfcc75805b45964
mod_leafnc_uncertainty <-summary(n1)$sigma
#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
lnpp_b <- summary(n1)$coefficients[1,1]/0.46 + summary(n1)$coefficients[2,1]*(vcmax25_df$vcmax25/LMA$myvar)
#therefore, the uncertainty of npp / gpp (defined as a here) = uncertainty b * (deriative a / deriative b)
mod_leafnc_uncertainty <- mod_leafnc_uncertainty *  exp(-lnpp_b) / ( (1 + exp(-lnpp_b)) ^2)

#NRE
mod_nre_uncertainty <- summary(nre_model)$sigma
nre_b <- summary(nre_model)$coefficients[1,1]  +
  summary(nre_model)$coefficients[2,1] * Tg$myvar + summary(nre_model)$coefficients[3,1] * log(vpd$myvar)
mod_nre_uncertainty <- mod_nre_uncertainty *  exp(-nre_b) / ( (1 + exp(-nre_b)) ^2)

#Uncertainty of GPP - assumed as standard deviation between observed and predicted gpp from FLUXNET (Stocker et al. 2020 GMD)
load("/Users/yunpeng/data/gpp_gmd/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata")
#Data in Euler is from: /cluster/work/climate/bestocke/data/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata 
#Data in my desktop is from: /Users/yunpeng/data/gpp_gmd/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata 
obs <-out_eval_FULL$gpp$fluxnet$data$meandf$obs
pred <- out_eval_FULL$gpp$fluxnet$data$meandf$mod
obs_pred <- as.data.frame(cbind(obs,pred))
obs_pred <- na.omit(obs_pred)
summary(obs_pred)
#analyse_modobs2(obs_pred,"pred","obs_pred", type = "points") # r2 is truly cloased to 0.7 - good 

# According to textbook:Uncertainty estimates obtained as standard deviations of repeated measurement results are called A type uncertainty estimates. 
# In this way,  “sample mean” corresponds to the true GPP, and the “observation” is the modelled GPP.
obs_pred$variance <- (obs_pred$obs - obs_pred$pred)^2
uncertainty_gpp <- sqrt(sum(obs_pred$variance)/nrow(obs_pred))
uncertainty_gpp
mean_gpp <- mean(obs_pred$pred,na.rm=TRUE)
mean_gpp
gpp_pft <- gpp_df$gpp
#gpp_pft[gpp_pft==0] <-NA
#summary(gpp_pft)

#now, calculate N uptake in the leaf = GPP * (ANPP/GPP) * (leafNPP/ANPP) * (leaf n/c) * (1-NRE)
uncertainty_lnf <- lnf_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                   (mod_anpp_uncertainty/(anpp_f/gpp_pft))^2 +
                                   (mod_lnpp_uncertainty/(lnpp_f/anpp_f))^2 +
                                   (mod_leafnc_uncertainty/leafcn_f)^2 +
                                   (mod_nre_uncertainty/(1-nre_f))^2)
summary(uncertainty_lnf/lnf_f)

#now, uncertainty of wood n uptake flux (assuming wood c/n is a constant = 100, without uncertainty) = GPP * (ANPP/GPP) * (1-leafNPP/ANPP) * (1/100)
uncertainty_wnf <- wnf_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (mod_anpp_uncertainty/(anpp_f/gpp_pft))^2 +
                                        (mod_lnpp_uncertainty/(1 - (lnpp_f/anpp_f)))^2)
summary(uncertainty_wnf/wnf_f)

#now, uncertainty of root n uptake flux (assuming root c/n is a constant = 94, without uncertainty)
#uncertainty of npp firstly
uncertainty_npp <- npp_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (mod_tnpp_uncertainty/(npp_f/gpp_pft))^2)
summary(uncertainty_npp/npp_f)

uncertainty_anpp <- anpp_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                           (mod_anpp_uncertainty/(anpp_f/gpp_pft))^2)
summary(uncertainty_anpp/anpp_f)

bnpp_gpp_uncertainty <- sqrt((mod_tnpp_uncertainty)^2 + (mod_anpp_uncertainty)^2)
uncertainty_bnpp <- bnpp_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                           (bnpp_gpp_uncertainty/(bnpp_f/gpp_pft))^2)
hist(uncertainty_bnpp/bnpp_f)
summary(uncertainty_bnpp/bnpp_f)

uncertainty_bnf <- (bnf_f/bnpp_f) * uncertainty_bnpp 
hist(uncertainty_bnf/bnf_f)
summary(uncertainty_bnf/bnf_f)

uncertainty_wnpp <- wnpp_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                           (mod_anpp_uncertainty/(anpp_f/gpp_pft))^2 +
                                           (mod_lnpp_uncertainty/(1 - (lnpp_f/anpp_f)))^2)
summary(uncertainty_wnpp/wnpp_f)

uncertainty_wnpp_gpp <- (sqrt((mod_anpp_uncertainty/(anpp_f/gpp_pft))^2 +
        (mod_lnpp_uncertainty/(1 - (lnpp_f/anpp_f)))^2))

summary(uncertainty_wnpp_gpp)

uncertainty_lnpp <- lnpp_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                           (mod_anpp_uncertainty/(anpp_f/gpp_pft))^2 +
                                           (mod_lnpp_uncertainty/(lnpp_f/anpp_f))^2)
summary(uncertainty_lnpp/lnpp_f)

uncertainty_nuptake <- sqrt(uncertainty_lnf^2 + uncertainty_wnf^2 + uncertainty_bnf^2)
summary(uncertainty_nuptake/nuptake_f)

uncertainty_forest <- as.data.frame(cbind(uncertainty_npp,uncertainty_anpp,uncertainty_bnpp,uncertainty_lnpp,uncertainty_wnpp,uncertainty_lnf,uncertainty_wnf,uncertainty_bnf,uncertainty_nuptake))
summary(uncertainty_forest)

#now, do the same for grassland
load(file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/tnpp_grass.RData")
load(file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/anpp_grass.RData")

#npp = 0.435 * gpp
uncertainty_grass_npp <- npp_g * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                                          (summary(tnpp_grass)$coef[1,2]/summary(tnpp_grass)$coef[1,1])^2)

uncertainty_grass_npp_gpp <- summary(tnpp_grass)$coef[1,2]

uncertainty_grass_anpp_gpp <- summary(anpp_grass)$coef[1,2]



#anpp = 0.228 * gpp
uncertainty_grass_anpp <- anpp_g * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                                             (summary(anpp_grass)$coef[1,2]/summary(anpp_grass)$coef[1,1])^2)

summary(uncertainty_grass_anpp/anpp_g)

#bnpp = gpp * (npp/gpp - anpp/gpp)
uncertainty_grass_bnpp <- sqrt((uncertainty_grass_npp)^2 + (uncertainty_grass_anpp)^2)
summary(uncertainty_grass_bnpp/bnpp_g)

#lnf = anpp * (1/18.0)* (1-NRE)
uncertainty_grass_lnf <- (lnf_g/anpp_g) * uncertainty_grass_anpp
summary(uncertainty_grass_lnf/lnf_g)

#bnf = bnpp * (1/41.0)
uncertainty_grass_bnf <- (bnf_g/bnpp_g) * uncertainty_grass_bnpp 
summary(uncertainty_grass_bnf/bnf_g)

uncertainty_grass_nuptake <- sqrt(uncertainty_grass_lnf^2 + uncertainty_grass_bnf^2)
summary(uncertainty_grass_nuptake/nuptake_g)

uncertainty_grass <- as.data.frame(cbind(uncertainty_grass_npp,uncertainty_grass_anpp,uncertainty_grass_bnpp,
                                         uncertainty_grass_lnf,uncertainty_grass_bnf,uncertainty_grass_nuptake))
summary(uncertainty_grass)


#uncertainty for forest
sum(uncertainty_npp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_anpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_bnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_lnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
#uncertainty of leaf npp/gpp

sum(mod_lnpp_uncertainty*(forest_percent *conversion)*available_grid2,na.rm=TRUE)

#uncertainty of wood npp/gpp
sum(uncertainty_wnpp_gpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)


sum(uncertainty_wnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_lnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_bnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_wnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_nuptake*(forest_percent *conversion)*available_grid2,na.rm=TRUE)

#grassland
sum(uncertainty_grass_npp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_grass_anpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_grass_bnpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_grass_lnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_grass_bnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)
sum(uncertainty_grass_nuptake*(grass_percent *conversion)*available_grid2,na.rm=TRUE)

#uncertainty of forest and grassland
uncertainty_gpp/mean_gpp
#uncertainty of npp/gpp
sqrt(sum(mod_tnpp_uncertainty*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + 
       sum(uncertainty_grass_npp_gpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
#uncertainty of npp
sqrt(sum(uncertainty_npp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_npp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
#uncertainty of anpp/gpp
sqrt(sum(mod_anpp_uncertainty*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + 
       sum(uncertainty_grass_anpp_gpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)

sqrt(sum(uncertainty_anpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_anpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_bnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_bnpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)

sqrt(sum(uncertainty_lnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_lnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_bnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_bnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_nuptake*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_nuptake*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)

#uncertainty of total nuptake and nue
aa <- sqrt((uncertainty_nuptake*forest_percent*available_grid2)^2 +
    (uncertainty_grass_nuptake*grass_percent*available_grid2)^2)

simulations$uncertainty_nuptake <- aa
simulations_predictors <- as.data.frame(cbind(all_predictors[,-c(9,10,11)],simulations))

simulations$uncertainty_nuptake_gpp <- simulations$uncertainty_nuptake/simulations$gpp

plot_map3(simulations[,c("lon","lat","uncertainty_nuptake")],
          varnam = "uncertainty_nuptake",
          latmin = -65, latmax = 85)
plot_map3(simulations[,c("lon","lat","uncertainty_nuptake_gpp")],
          varnam = "uncertainty_nuptake_gpp",
          latmin = -65, latmax = 85)


'''
######plot all maps######
#run a loop function to output map - and output total values
for (i in 4:ncol(all_maps)){
  varname <- names(all_maps)[i]
  total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
  plot_map3(all_maps[,c("lon","lat",varname)],
            varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
            latmin = -65, latmax = 85)
  ggsave(paste("/Users/yunpeng/data/output/output_map/no_points/",varname,".jpg",sep=""))
}

#now, read all_points from forest, grass, leafnc, nre
#1. gpp
NPP_forest <- read.csv("/Users/yunpeng/data/NPP_final/NPP_validation.csv")
NPP_grassland <- read.csv("/Users/yunpeng/data/NPP_Grassland_final/NPP_grass_validation.csv")

gpp_f <- (NPP_forest %>% filter(GPP>0) %>% filter(pred_gpp_c3>0))[,c("lon","lat")]
gpp_g <- (NPP_grassland %>% filter(GPP>0) %>% filter(weightedgpp_measured_c3>0))[,c("lon","lat")]

i=4
varname <- names(all_maps)[i]
varname
total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=gpp_f,aes(lon,lat),col="red")+
  geom_point(data=gpp_g,aes(lon,lat),col="blue")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))

#2.tnpp
npp_f <- (NPP_forest %>% filter(TNPP_1>0) %>% filter(pred_npp>0))[,c("lon","lat")]
npp_g <- (NPP_grassland %>% filter(TNPP_1>0) %>% filter(pred_npp>0))[,c("lon","lat")]

i=5
varname <- names(all_maps)[i]
varname
total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
          varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
          latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=npp_f,aes(lon,lat),col="red")+
  geom_point(data=npp_g,aes(lon,lat),col="blue")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))

#3.anpp
anpp_f <- (NPP_forest %>% filter(ANPP_2>0) %>% filter(pred_anpp>0))[,c("lon","lat")]
anpp_g <- (NPP_grassland %>% filter(ANPP_2>0) %>% filter(pred_anpp>0))[,c("lon","lat")]

i=8
varname <- names(all_maps)[i]
varname
total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=anpp_f,aes(lon,lat),col="red")+
  geom_point(data=anpp_g,aes(lon,lat),col="blue")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))

#4.leaf npp
lnpp_f <- (NPP_forest %>% filter(NPP.foliage>0) %>% filter(pred_lnpp>0))[,c("lon","lat")]

i=14
varname <- names(all_maps)[i]
varname
total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=lnpp_f,aes(lon,lat),col="red")  
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))

# 5. bnpp
bnpp_f <- (NPP_forest %>% filter(BNPP_1>0) %>% filter(pred_bnpp>0))[,c("lon","lat")]
bnpp_g <- (NPP_grassland %>% filter(BNPP_1>0) %>% filter(pred_bnpp>0))[,c("lon","lat")]

i=11
varname <- names(all_maps)[i]
varname
total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=bnpp_f,aes(lon,lat),col="red")+
  geom_point(data=bnpp_g,aes(lon,lat),col="blue")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))

# 6. laefnc
SP_input <- read.csv(file="/Users/yunpeng/data/leaf_traits/combined_leaf_traits.csv") #new one 
SP_input2 <- SP_input[,c("lat","lon","z","Vcmax25","narea","lma")]
sitemean <- aggregate(SP_input2,by=list(SP_input2$lon,SP_input2$lat), FUN=mean, na.rm=TRUE) 
sitemean$pred_leafn <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax25/sitemean$lma
sitemean$obs_leafn <- sitemean$narea/sitemean$lma

laefnc_f <- (sitemean %>% filter(pred_leafn>0) %>% filter(obs_leafn>0))[,c("lon","lat")]

i=16
varname <- names(all_maps)[i]
varname
total_value <- round(mean(all_maps[,i],na.rm=TRUE),3)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste("leaf N/C", ":", total_value, "(assuming C% as global constant = 48%)", sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=laefnc_f,aes(lon,lat),col="red")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))

# 7. NRE
NRE_validation <- read.csv(file="/Users/yunpeng/data/NPP_final/NRE_validation.csv") #new one 

nre_f <- (NRE_validation %>% filter(pred_nre>0) %>% filter(NRE>0))[,c("lon","lat")]

i=17
varname <- names(all_maps)[i]
varname
total_value <- round(mean(all_maps[,i],na.rm=TRUE),3)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste("N resorption efficiency", ":", total_value, sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=nre_f,aes(lon,lat),col="red")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))


# 8. lnf
lnf_f <- (NPP_forest %>% filter(pred_lnf>0) %>% filter(lnf_obs_org>0))[,c("lon","lat")]
lnf_g <- (NPP_grassland %>% filter(pred_lnf>0) %>% filter(lnf_obs_org>0))[,c("lon","lat")]

i=19
varname <- names(all_maps)[i]
varname
total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=lnf_f,aes(lon,lat),col="red")+
  geom_point(data=lnf_g,aes(lon,lat),col="blue")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))


# 9. nuptake
Nmin_validation <- read.csv("/Users/yunpeng/data/NPP_final/Nmin_validation.csv")
nuptake_f <- (Nmin_validation %>% filter(pred_nuptake>0) %>% filter(obs_nuptake>0))[,c("lon","lat")]

i=25
varname <- names(all_maps)[i]
varname
total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
gg <- plot_map3(all_maps[,c("lon","lat",varname)],
                varnam = varname,plot_title = paste(varname, ":", total_value, "Pg/yr", sep=" " ),
                latmin = -65, latmax = 85,combine=FALSE)
gg$ggmap +
  geom_point(data=nuptake_f,aes(lon,lat),col="red")
ggsave(paste("/Users/yunpeng/data/output/output_map/",varname,".jpg",sep=""))


