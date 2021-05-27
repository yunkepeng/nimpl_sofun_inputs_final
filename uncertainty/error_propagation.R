rm(list=ls())
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

#### Forest
#1. In our path (with multiple years data), identify which is the first year and end year of those files
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

npp_df <- inputnc("npp",1982,2011)

anpp_df <- inputnc("anpp",1982,2011)

bnpp_df <- inputnc("bnpp",1982,2011)

lnpp_df <- inputnc("lnpp",1982,2011)

wnpp_df <- inputnc("wnpp",1982,2011)

leafcn_df <- inputnc("leafcn",1982,2011) # this is actually leaf n/c. 

lnf_df <- inputnc("lnf",1982,2011) 

wnf_df <- inputnc("wnf",1982,2011) 

bnf_df <- inputnc("bnf",1982,2011) 

nuptake_df <- inputnc("nuptake",1982,2011) 

nre_df <- inputnc("nre",1982,2011)

#### Grassland
##In our path (with multiple years data), identify which is the first year and end year of those files
firstyr_data <- 1982 # In data file, which is the first year
endyr_data <- 2011 # In data file, which is the last year
location <- "/Users/yunpeng/data/output/latest_grass/"
alloutput_list <- list.files(location,full.names = T)

grass_inputnc <- function(name,start_year,end_year){
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
    if (name == "npp_grass"){
      nc <- read_nc_onefile(alloutput_list[grepl("a.npp_grass.nc", list.files(location,full.names = T))][i-firstyr_data+1]) #we only rely this to filter npp.nc file...
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

npp_grass_df <- grass_inputnc("npp_grass",1982,2011)

anpp_grass_df <- grass_inputnc("anpp_grass",1982,2011)

bnpp_grass_df <- grass_inputnc("bnpp_grass",1982,2011)

lnf_grass_df <- grass_inputnc("lnf_grass",1982,2011) 

bnf_grass_df <- grass_inputnc("bnf_grass",1982,2011) 

nuptake_grass_df <- grass_inputnc("nuptake_grass",1982,2011) 

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

#save.image(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/uncertainty/uncertainty.Rdata")

rm(list=ls())
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/uncertainty/uncertainty.Rdata")

#for consistency of gpp, npp....nuptake in forest and grassland check, between FORTRAN and R: see consistency_check.R - as located in the same file
#after check in the code - go ahead
#Uncertainty of GPP - assumed as standard deviation between observed and predicted gpp from FLUXNET (Stocker et al. 2020 GMD)
load("/Users/yunpeng/data/gpp_gmd/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata")
#Data in Euler is from: /cluster/work/climate/bestocke/data/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata 
#Data in my desktop is from: /Users/yunpeng/data/gpp_gmd/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata 
obs <-out_eval_FULL$gpp$fluxnet$data$meandf$obs
pred <- out_eval_FULL$gpp$fluxnet$data$meandf$mod
obs_pred <- as.data.frame(cbind(obs,pred))
obs_pred <- na.omit(obs_pred)
summary(obs_pred)
analyse_modobs2(obs_pred,"pred","obs_pred", type = "points") # r2 is truly cloased to 0.7 - good 

# According to textbook:Uncertainty estimates obtained as standard deviations of repeated measurement results are called A type uncertainty estimates. 
# In this way,  “sample mean” corresponds to the true GPP, and the “observation” is the modelled GPP.
obs_pred$variance <- (obs_pred$obs - obs_pred$pred)^2
uncertainty_gpp <- sqrt(sum(obs_pred$variance)/nrow(obs_pred))
uncertainty_gpp
mean_gpp <- mean(obs_pred$pred,na.rm=TRUE)

#Error propagation of NPP - firstly - load coefficients from all models
#firstly, load all forest models
load("/Users/yunpeng/data/NPP_final/statistical_model/mod_tnpp.RData")
summary(mod_tnpp)
load("/Users/yunpeng/data/NPP_final/statistical_model/mod_anpp.RData")
summary(mod_anpp)
load("/Users/yunpeng/data/NPP_final/statistical_model/mod_lnpp.RData")
summary(mod_lnpp)
load("/Users/yunpeng/data/NPP_final/statistical_model/nmass.RData")
summary(n1)
load("/Users/yunpeng/data/NPP_final/statistical_model/nre_model.RData")
summary(nre_model)

## TNPP/GPP
#using standard error method to calculate uncertainty of whole regression
# for all logit function model (mod_tnpp, mod_anpp, mod_lnpp, nre_model): a = 1/(1+exp(-b)), where a is the ratio (e.g. npp/gpp) and b is the regression (e.g. mod =tnpp)
#we calculate uncertainty of b firstly, based on each regression
# npp/gpp uncertainty
mod_tnpp_uncertainty <- summary(mod_tnpp)$sigma  #random factor: a population parameter - estimated automatically
mod_tnpp_uncertainty
#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
tnpp_b <- summary(mod_tnpp)$coefficients[1,1] + summary(mod_tnpp)$coefficients[2,1] * log(CNrt[,3])  +
  summary(mod_tnpp)$coefficients[3,1] * log(age[,3]) + summary(mod_tnpp)$coefficients[4,1] * fAPAR[,3]

#therefore, the uncertainty of npp / gpp (defined as a here) = uncertainty b * (deriative a / deriative b)
mod_tnpp_uncertainty <- mod_tnpp_uncertainty *  exp(-tnpp_b) / ( (1 + exp(-tnpp_b)) ^2)

## ANPP/GPP
mod_anpp_uncertainty <- summary(mod_anpp)$sigma

#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
anpp_b <- summary(mod_anpp)$coefficients[1,1] + summary(mod_anpp)$coefficients[2,1] * log(CNrt[,3])  +
  summary(mod_anpp)$coefficients[3,1] * log(age[,3]) + summary(mod_anpp)$coefficients[4,1] * fAPAR[,3]

#therefore, the uncertainty of npp / gpp (defined as a here) = uncertainty b * (deriative a / deriative b)
mod_anpp_uncertainty <- mod_anpp_uncertainty *  exp(-anpp_b) / ( (1 + exp(-anpp_b)) ^2)

## LNPP/ANPP
mod_lnpp_uncertainty <- summary(mod_lnpp)$sigma

#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
lnpp_b <- summary(mod_lnpp)$coefficients[1,1] + summary(mod_lnpp)$coefficients[2,1] * log(PPFD[,3])  +
  summary(mod_lnpp)$coefficients[3,1] * Tg[,3] + summary(mod_lnpp)$coefficients[4,1] * log(vpd[,3])

#therefore, the uncertainty of npp / gpp (defined as a here) = uncertainty b * (deriative a / deriative b)
mod_lnpp_uncertainty <- mod_lnpp_uncertainty *  exp(-lnpp_b) / ( (1 + exp(-lnpp_b)) ^2)

#leaf n/c - ALREADY checked that inputted vcmax25_df + LMA calculated from (1) fortran and (2) R for predicting leaf n/c, could output the same prediction for leaf n/c Please note! using a.vcmax25.nc rather than annualvcmax25.nc. See difference in: https://www.notion.so/computationales/annualvcmax25-vs-a-vcmax25-282bddd63bba4d7b9cfcc75805b45964
mod_leafnc_uncertainty <-summary(n1)$sigma

#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
lnpp_b <- summary(n1)$coefficients[1,1]/0.46 + summary(n1)$coefficients[2,1]*(vcmax25_df[,4]/LMA[,3])

#therefore, the uncertainty of npp / gpp (defined as a here) = uncertainty b * (deriative a / deriative b)
mod_leafnc_uncertainty <- mod_leafnc_uncertainty *  exp(-lnpp_b) / ( (1 + exp(-lnpp_b)) ^2)

#NRE
mod_nre_uncertainty <- summary(nre_model)$sigma

#now we need to calculate deriative a / deriative b, based on a = 1/(1+exp(-b)). After calculation it is: exp(-b) / ( (1 + exp(-b)) ^2)
lnpp_b <- summary(nre_model)$coefficients[1,1]  +
  summary(nre_model)$coefficients[2,1] * Tg[,3] + summary(nre_model)$coefficients[3,1] * log(vpd[,3])

#therefore, the uncertainty of npp / gpp (defined as a here) = uncertainty b * (deriative a / deriative b)
mod_nre_uncertainty <- mod_nre_uncertainty *  exp(-lnpp_b) / ( (1 + exp(-lnpp_b)) ^2)

#now, calculate N uptake in the leaf = GPP * (ANPP/GPP) * (leafNPP/ANPP) * (leaf n/c) * (1-NRE)
uncertainty_lnf <- lnf_df$lnf * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (mod_anpp_uncertainty/(anpp_df$anpp/gpp_df$gpp))^2 +
                                        (mod_lnpp_uncertainty/(lnpp_df$lnpp/anpp_df$anpp))^2 +
                                        (mod_leafnc_uncertainty/leafcn_df$leafcn)^2 +
                                        (mod_nre_uncertainty/(1-nre_df$nre))^2)
summary(uncertainty_lnf)
summary(uncertainty_lnf/lnf_df$lnf)

#now, uncertainty of wood n uptake flux (assuming wood c/n is a constant = 100, without uncertainty) = GPP * (ANPP/GPP) * (1-leafNPP/ANPP) * (1/100)
uncertainty_wnf <- wnf_df$wnf * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (mod_anpp_uncertainty/(anpp_df$anpp/gpp_df$gpp))^2 +
                                        (mod_lnpp_uncertainty/(1 - (lnpp_df$lnpp/anpp_df$anpp)))^2)
summary(uncertainty_wnf/wnf_df$wnf)

#now, uncertainty of root n uptake flux (assuming root c/n is a constant = 94, without uncertainty)
#uncertainty of npp firstly
uncertainty_npp <- npp_df$npp * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (mod_tnpp_uncertainty/(npp_df$npp/gpp_df$gpp))^2)
summary(uncertainty_npp/npp_df$npp)

uncertainty_anpp <- anpp_df$anpp * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                           (mod_anpp_uncertainty/(anpp_df$anpp/gpp_df$gpp))^2)
summary(uncertainty_anpp/anpp_df$anpp)

bnpp_gpp_uncertainty <- sqrt((mod_tnpp_uncertainty)^2 + (mod_anpp_uncertainty)^2)
uncertainty_bnpp <- bnpp_df$bnpp * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (bnpp_gpp_uncertainty/(bnpp_df$bnpp/gpp_df$gpp))^2)
hist(uncertainty_bnpp/bnpp_df$bnpp)
summary(uncertainty_bnpp/bnpp_df$bnpp)

uncertainty_bnf <- (bnf_df$bnf/bnpp_df$bnpp) * uncertainty_bnpp 
summary(uncertainty_bnf/bnf_df$bnf)

uncertainty_wnpp <- wnpp_df$wnpp * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                           (mod_anpp_uncertainty/(anpp_df$anpp/gpp_df$gpp))^2 +
                                           (mod_lnpp_uncertainty/(1 - (lnpp_df$lnpp/anpp_df$anpp)))^2)
summary(uncertainty_wnpp/wnpp_df$wnpp)

uncertainty_lnpp <- lnpp_df$lnpp * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (mod_anpp_uncertainty/(anpp_df$anpp/gpp_df$gpp))^2 +
                                        (mod_lnpp_uncertainty/(lnpp_df$lnpp/anpp_df$anpp))^2)
summary(uncertainty_lnpp/lnpp_df$lnpp)

uncertainty_nuptake <- sqrt(uncertainty_lnf^2 + uncertainty_wnf^2 + uncertainty_bnf^2)
summary(uncertainty_nuptake/nuptake_df$nuptake)

uncertainty_forest <- as.data.frame(cbind(uncertainty_npp,uncertainty_anpp,uncertainty_bnpp,uncertainty_lnpp,uncertainty_wnpp,uncertainty_lnf,uncertainty_wnf,uncertainty_bnf,uncertainty_nuptake))
summary(uncertainty_forest)

#now, do the same for grassland
load(file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/tnpp_grass.RData")
load(file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/anpp_grass.RData")

#npp = 0.435 * gpp
uncertainty_grass_npp <- npp_grass_df$npp_grass * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                        (summary(tnpp_grass)$coef[1,2]/summary(tnpp_grass)$coef[1,1])^2)

summary(uncertainty_grass_npp/npp_grass_df$npp_grass)

#anpp = 0.228 * gpp
uncertainty_grass_anpp <- anpp_grass_df$anpp_grass * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                                          (summary(anpp_grass)$coef[1,2]/summary(anpp_grass)$coef[1,1])^2)

summary(uncertainty_grass_anpp/anpp_grass_df$anpp_grass)

#bnpp = gpp * (npp/gpp - anpp/gpp)
uncertainty_grass_bnpp <- sqrt((uncertainty_grass_npp)^2 + (uncertainty_grass_anpp)^2)
summary(uncertainty_grass_bnpp/bnpp_grass_df$bnpp_grass)

#lnf = anpp * (1/18.0)* (1-NRE)
uncertainty_grass_lnf <- lnf_grass_df$lnf_grass * sqrt( (uncertainty_grass_anpp/anpp_grass_df$anpp_grass)^2 +
                                        (mod_nre_uncertainty/(1-nre_df$nre))^2)
summary(uncertainty_grass_lnf/lnf_grass_df$lnf_grass)

#bnf = bnpp * (1/41.0)
uncertainty_grass_bnf <- (bnf_grass_df$bnf_grass/bnpp_grass_df$bnpp_grass) * uncertainty_grass_bnpp 
summary(uncertainty_grass_bnf/bnf_grass_df$bnf_grass)

uncertainty_grass_nuptake <- sqrt(uncertainty_grass_lnf^2 + uncertainty_grass_bnf^2)
summary(uncertainty_grass_nuptake/nuptake_grass_df$nuptake_grass)

uncertainty_grass <- as.data.frame(cbind(uncertainty_grass_npp,uncertainty_grass_anpp,uncertainty_grass_bnpp,
                                         uncertainty_grass_lnf,uncertainty_grass_bnf,uncertainty_grass_nuptake))
summary(uncertainty_grass)

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

#check thousands of points with NA
lonlat <- gpp_df[,c("lon","lat")]

####now, firstly calculate gridded map sum

#area_m2 to show each grid's area in m2
calc_area <- function( lat, dx=1, dy=1 ){
  r_earth <- 6370499.317638  # to be consistent with how Ferret calculates areas of spheres (https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2016/msg00155.html)
  area <- 4 * r_earth^2 * 0.5 * dx * pi/180 * cos( abs(lat) * pi/180 ) * sin( 0.5 * dy * pi/180 )
  return(area)
}
area_m2 <- calc_area(lonlat$lat,0.5,0.5)

#fland - to show each grid's land cover percentage
nc <- read_nc_onefile("/Users/yunpeng/data/fland/global.fland.nc") #Input nc
output_fland <- nc_to_df(nc, varnam = "fland")
fland <- output_fland$myvar

#also, remove the grids with fortran outputted value as na 
all_predictors2 <- as.data.frame(cbind(lonlat,npp_df[,4],lnpp_df[,4],lnf_df[,4],gpp_df[,4]))
all_predictors2$avail <- 1
names(all_predictors2) <- c("lon","lat","npp","lnpp","lnf","gpp","avail")
#convert all 0 values to NA (they are either due to blank values for predictors, or due to gpp =0)
all_predictors2$avail[(all_predictors2$npp)==0 | (all_predictors2$lnpp)==0| (all_predictors2$lnf)==0] <- NA
available_grids2 <- all_predictors2$avail

dim(subset(all_predictors2,is.na(avail)==TRUE)) # 12720,but actually 5245 points were NA - other 7475 points were gpp = 0

library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
#red points were prediction fields =NA
points(subset(all_predictors2,is.na(avail)==TRUE)$lon,subset(all_predictors2,is.na(avail)==TRUE)$lat, col="red", pch=16,cex=1)

#blue points were gpp = 0
points(subset(gpp_df,gpp==0)$lon,subset(gpp_df,gpp==0)$lat, col="blue", pch=16,cex=1)
#we keep blue points still equal to 0, but red points converted to 0

#now, examine some strange points that enttered in Fortran as NA (0), but in R it works
all_predictors <- as.data.frame(cbind(lonlat,Tg$myvar,PPFD$myvar,vpd$myvar,
                                      alpha$myvar,fAPAR$myvar,age$myvar,
                                      CNrt$myvar,LMA$myvar))
all_predictors$avail_fortran <-available_grids2
all_predictors$gpp <- gpp_df$gpp
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
#!!!!!!red points were prediction fields =NA but in R it works
summary(subset(all_predictors,is.na(avail_fortran)==TRUE & gpp>0)) #see..not because of prediction fields
points(subset(all_predictors,is.na(avail_fortran)==TRUE & gpp>0)$lon,subset(all_predictors,is.na(avail_fortran)==TRUE & gpp>0)$lat, col="red", pch=16,cex=1)



#now, input all files - weighted sum - and also multiply with available grids to filter NA grids in all map (except for gpp)
gpp_final <- gpp_df[,4]# forest only

npp_final_forest <- available_grids2 * (npp_df[,4] * forest_percent)
npp_final_grass <- available_grids2 * (npp_grass_df[,4] * grass_percent)
npp_final <- available_grids2 * (npp_df[,4] * forest_percent + npp_grass_df[,4] * grass_percent)

anpp_final_forest <- available_grids2 *  (anpp_df[,4] * forest_percent)
anpp_final_grass <- available_grids2 *  (anpp_grass_df[,4] * grass_percent)
anpp_final <- available_grids2 *  (anpp_df[,4] * forest_percent + anpp_grass_df[,4] * grass_percent)

bnpp_final_forest <- available_grids2 *  (bnpp_df[,4] * forest_percent)
bnpp_final_grass <- available_grids2 *  (bnpp_grass_df[,4] * grass_percent)
bnpp_final <- available_grids2 *  (bnpp_df[,4] * forest_percent + bnpp_grass_df[,4] * grass_percent)

lnpp_final_forest <- available_grids2 * lnpp_df[,4] * forest_percent  #forest only

wnpp_final_forest <- available_grids2 *  wnpp_df[,4] * forest_percent  #forest only

lnf_final_forest <- available_grids2 *  (lnf_df[,4] * forest_percent)
lnf_final_grass <- available_grids2 *  (lnf_grass_df[,4] * grass_percent)
lnf_final <- available_grids2 *  (lnf_df[,4] * forest_percent + lnf_grass_df[,4] * grass_percent)

bnf_final_forest <- available_grids2 *  (bnf_df[,4] * forest_percent)
bnf_final_grass <- available_grids2 *  (bnf_grass_df[,4] * grass_percent)
bnf_final <- available_grids2 * (bnf_df[,4] * forest_percent + bnf_grass_df[,4] * grass_percent)

wnf_final_forest <- available_grids2 *wnf_df[,4] * forest_percent  #forest only

nuptake_final_forest <- available_grids2 *  (nuptake_df[,4] * forest_percent)
nuptake_final_grass <- available_grids2 *  (nuptake_grass_df[,4] * grass_percent)
nuptake_final <- available_grids2 * (nuptake_df[,4] * forest_percent + nuptake_grass_df[,4] * grass_percent)

final_all <- as.data.frame(cbind(lonlat,gpp_final,npp_final_forest,npp_final_grass,npp_final,anpp_final_forest,anpp_final_grass,anpp_final,
                                 bnpp_final_forest,bnpp_final_grass,bnpp_final,lnpp_final_forest,wnpp_final_forest,lnf_final_forest,lnf_final_grass,lnf_final,
                                 bnf_final_forest,bnf_final_grass,bnf_final,wnf_final_forest,nuptake_final_forest,nuptake_final_grass,nuptake_final))
dim(final_all)
head(final_all)



load(file = "~/yunkepeng/nimpl_sofun_inputs/forest/New_Nuptake_site_simulation.Rdata")
gg <- plot_map3(final_all[,c("lon","lat","nuptake_final")], 
                varnam = "nuptake_final",plot_title = "Total N uptake (gN/m2/yr)",
                latmin = -65, latmax = 85, combine = FALSE)

gg$ggmap + geom_point(data=subset(Nmin_all,pred_nuptake>0),aes(lon,lat),size=3,col="red")
gg$gglegend 

#calculate outputted values
conversion <- area_m2 * fland /1e+15 #this is conversion factor --> (1) inclulding grid's area and land cover and (2) converting from gC/m2/yr to PgC/yr - just add sum then it works

#lastly, convert all gpp = 0's grid - for their npp also = 0, not na
for (i in 3:(ncol(final_all))){
  final_all[,i][final_all$gpp_final==0] <- 0
  print(names(final_all)[i])
  print( sum(final_all[,i]*conversion,na.rm=TRUE))
}

#uncertainty of gpp

dim(subset(final_all,gpp_final==0))
dim(subset(final_all,npp_final==0))
dim(subset(final_all,anpp_final==0))
dim(subset(final_all,nuptake_final==0))


sum(uncertainty_gpp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
#uncertainty for forest
sum(uncertainty_npp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_anpp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_bnpp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_lnpp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_wnpp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_lnf*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_bnf*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_wnf*(forest_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_nuptake*(forest_percent *conversion)*available_grids2,na.rm=TRUE)

#grassland
sum(uncertainty_grass_npp*(grass_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_grass_anpp*(grass_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_grass_bnpp*(grass_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_grass_lnf*(grass_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_grass_bnf*(grass_percent *conversion)*available_grids2,na.rm=TRUE)
sum(uncertainty_grass_nuptake*(grass_percent *conversion)*available_grids2,na.rm=TRUE)

#uncertainty of forest and grassland
sqrt(sum(uncertainty_npp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)^2 + sum(uncertainty_grass_npp*(grass_percent *conversion)*available_grids2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_anpp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)^2 + sum(uncertainty_grass_anpp*(grass_percent *conversion)*available_grids2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_bnpp*(forest_percent *conversion)*available_grids2,na.rm=TRUE)^2 + sum(uncertainty_grass_bnpp*(grass_percent *conversion)*available_grids2,na.rm=TRUE)^2)

sqrt(sum(uncertainty_lnf*(forest_percent *conversion)*available_grids2,na.rm=TRUE)^2 + sum(uncertainty_grass_lnf*(grass_percent *conversion)*available_grids2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_bnf*(forest_percent *conversion)*available_grids2,na.rm=TRUE)^2 + sum(uncertainty_grass_bnf*(grass_percent *conversion)*available_grids2,na.rm=TRUE)^2)
sqrt(sum(uncertainty_nuptake*(forest_percent *conversion)*available_grids2,na.rm=TRUE)^2 + sum(uncertainty_grass_nuptake*(grass_percent *conversion)*available_grids2,na.rm=TRUE)^2)



#now, do the primary controls of N uptake
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

#CNrt$myvar <- mean(CNrt$myvar,na.rm=TRUE)
#age$myvar <- mean(age$myvar,na.rm=TRUE)
#fAPAR$myvar <- mean(fAPAR$myvar,na.rm=TRUE)
#PPFD$myvar <- mean(PPFD$myvar,na.rm=TRUE)
#vpd$myvar <- mean(vpd$myvar,na.rm=TRUE)
#Tg$myvar <- mean(Tg$myvar,na.rm=TRUE)
#vcmax25_df$vcmax25 <- 118.5
#LMA$myvar <- mean(LMA$myvar,na.rm=TRUE)

npp_f <- gpp_df$gpp * (1/(1 + exp(-(-0.36075 * log(CNrt$myvar) -0.16213 * log(age$myvar) + 
                                             0.72793 * fAPAR$myvar + 0.57014))))

anpp_f <- gpp_df$gpp * (1/(1 + exp(-(-0.55151 * log(CNrt$myvar) -0.20050 * log(age$myvar) + 
                                               1.06611 * fAPAR$myvar+ 0.35817))))

bnpp_f <- npp_f-anpp_f

lnpp_f <- anpp_df$anpp * (1/(1 + exp(-(0.97093* log(PPFD$myvar) +
                                                  0.06453 * (Tg$myvar) 
                                                -0.80397 * log(vpd$myvar)
                                                -7.47165))))

wnpp_f <- anpp_f - lnpp_f

leafcn_f <- (0.01599/0.46) + (0.005992/0.46) *vcmax25_df$vcmax25/LMA$myvar

nre_f <- (1/(1+exp(-(-0.064460 *Tg$myvar + 0.402850 * log(vpd$myvar) + 1.368935))))

lnf_f <- (1-nre_f)* leafcn_f * lnpp_f

wnf_f <- wnpp_f/100

bnf_f <- bnpp_f/94

nuptake_f <- lnf_f + wnf_f + bnf_f

###grassland
npp_g <- gpp_df$gpp *0.435
anpp_g <- gpp_df$gpp *0.228
bnpp_g <- npp_g - anpp_g
lnf_g <-  anpp_g *(1/18)*(1-(1/(1+exp(-(-0.064460 *Tg$myvar + 0.402850 * log(vpd$myvar) + 1.368935)))))
bnf_g <- bnpp_g *(1/41)
nuptake_g <- lnf_g + bnf_g

all_predictors <- as.data.frame(cbind(Tg$myvar,PPFD$myvar,vpd$myvar,
                                      alpha$myvar,fAPAR$myvar,age$myvar,
                                      CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25))

final_r <- as.data.frame(cbind(gpp_df,nuptake_f,nuptake_g,all_predictors))

final_r$total_uptake <- (forest_percent*final_r$nuptake_f + grass_percent*final_r$nuptake_g)

names(final_r) <- c("lon","lat","z","gpp","nuptake_f","nuptake_g","Tg","PPFD","vpd","alpha","fAPAR","age","CNrt","LMA","vcmax25","total_uptake")

a1 <- (lm(total_uptake~Tg+PPFD+vpd+fAPAR+CNrt+age+LMA+vcmax25,data=final_r))
library(visreg)
#visreg(a1)

anova(a1)
af <- anova(a1)
afss <- af$"Sum Sq"
options(scipen = 999)
print(cbind(af,PctExp=afss/sum(afss)*100))
