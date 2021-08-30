#uncertainty and estimated values
rm(list=ls())
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
library(rbeni)
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)

###1. load all prediction fields map (details in ~/yunkepeng/nimpl_sofun_inputs_final/Prediction_field), gpp and vcmax25 (derived from SOFUN/yunkebranch).
firstyr_data <- 1982 # In data file, which is the first year
endyr_data <- 2011 # In data file, which is the last year
location <- "~/data/output/latest_forest/"
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
#alpha not used...
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
nre_constant_grass <- 0.69

#input all regressions 

#firstly, load all forest models
load("~/data/NPP_final/statistical_model/mod_tnpp.RData")
summary(mod_tnpp)
load("~/data/NPP_final/statistical_model/mod_anpp.RData")
summary(mod_anpp)
load("~/data/NPP_final/statistical_model/mod_lnpp.RData")
summary(mod_lnpp)
load("~/data/NPP_final/statistical_model/nmass.RData")
summary(n1)
load("~/data/NPP_final/statistical_model/nre_model_forest.RData")
summary(nre_model)

#now, do the same for grassland
load(file = "~/data/NPP_grassland_final/statistical_model/tnpp_grass.RData")
mod_tnpp_grass<- tnpp_grass
summary(mod_tnpp_grass)
load(file = "~/data/NPP_grassland_final/statistical_model/anpp_grass.RData")
mod_anpp_grass <- anpp_grass
summary(mod_anpp_grass)

###2. run global simulations
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

leafnc_f <- (summary(n1)$coef[1,1]/0.46) + 
  (summary(n1)$coef[2,1]/0.46) *vcmax25_df$vcmax25/LMA$myvar
#0.46 is constant Cmass

nre_f <- (1/(1+exp(-(summary(nre_model)$coef[1,1]+
                       summary(nre_model)$coef[2,1] *Tg$myvar + 
                       summary(nre_model)$coef[3,1] * log(vpd$myvar)))))

lnf_f <- (1-nre_f)* leafnc_f * lnpp_f

wnf_f <- wnpp_f/100
#100 is constant wood c/n

bnf_f <- bnpp_f/94
#94 is constant root c/n

nuptake_f <- lnf_f + wnf_f + bnf_f

#grass
npp_g <- gpp_df$gpp * summary(mod_tnpp_grass)$coef[1,1]
anpp_g <- gpp_df$gpp * summary(mod_anpp_grass)$coef[1,1]
bnpp_g <- npp_g-anpp_g

leafnc_g <- 1/18

nre_g <- 0.69

lnf_g <- anpp_g*leafnc_g*(1-nre_g)

bnf_g <- bnpp_g *(1/41)
#41 is constant root c/n
nuptake_g <- lnf_g + bnf_g

###2. input land cover
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


#final calculation - now divide into forest, grassland and pft
#available_grid2 here was used as a list of data to identify if a grid is available (=1) or any prediction fields shown as NA 
npp_pft <- available_grid2* (npp_f*forest_percent +npp_g*grass_percent)
npp_forest <- available_grid2* (npp_f*forest_percent)
npp_grass <- available_grid2* (npp_g*grass_percent)

anpp_pft <- available_grid2*(anpp_f*forest_percent +anpp_g*grass_percent)
anpp_forest <- available_grid2* (anpp_f*forest_percent)
anpp_grass <- available_grid2* (anpp_g*grass_percent)

lnpp_forest <- available_grid2*lnpp_f*forest_percent

wnpp_forest <- available_grid2*wnpp_f*forest_percent

bnpp_pft <- available_grid2*(bnpp_f*forest_percent +bnpp_g*grass_percent)
bnpp_forest <- available_grid2* (bnpp_f*forest_percent)
bnpp_grass <- available_grid2* (bnpp_g*grass_percent)

leafcn_pft <- 1/(available_grid2*(leafnc_f*forest_percent +leafnc_g*grass_percent))
leafcn_forest <- 1/available_grid2*leafnc_f*forest_percent
leafcn_grassland <- 1/available_grid2*leafnc_g*grass_percent
summary(leafcn_pft)

nre_pft <- available_grid2*(nre_f*forest_percent +nre_g*grass_percent)
nre_forest <-  available_grid2*nre_f*forest_percent
nre_grassland <- available_grid2*nre_g*grass_percent
summary(nre_pft)

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
                                anpp_pft,anpp_forest,anpp_grass,
                                bnpp_pft,bnpp_forest,bnpp_grass,
                                lnpp_forest,wnpp_forest,wnf_forest,
                                leafcn_pft,leafcn_forest,leafcn_grassland,
                                nre_pft,nre_forest,nre_grassland,
                                lnf_pft,lnf_forest,lnf_grass,
                                bnf_pft,bnf_forest,bnf_grass,
                                nuptake_pft,nuptake_forest,nuptake_grass))

summary(all_maps)

#####area_m2 to show each grid's area in m2
calc_area <- function( lat, dx=1, dy=1 ){
  r_earth <- 6370499.317638  # to be consistent with how Ferret calculates areas of spheres (https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2016/msg00155.html)
  area <- 4 * r_earth^2 * 0.5 * dx * pi/180 * cos( abs(lat) * pi/180 ) * sin( 0.5 * dy * pi/180 )
  return(area)
}
lonlat <- gpp_df[,c("lon","lat")]
area_m2 <- calc_area(lonlat$lat,0.5,0.5)
#fland - to show each grid's land cover percentage
nc <- read_nc_onefile("~/data/fland/global.fland.nc") #Input nc
output_fland <- nc_to_df(nc, varnam = "fland")
fland <- output_fland$myvar
#include conversion factor (from g to Pg)
conversion <- area_m2 * fland /1e+15

#####Table S1
#just print estimations of all values
for (i in 4:31){
  varname <- names(all_maps)[i]
  total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
  print(varname)
  print(total_value)
}


#Uncertainty of GPP - assumed as standard deviation between observed and predicted gpp from FLUXNET (Stocker et al. 2020 GMD)
load("~/data/gpp_gmd/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata")
#Data in Euler is from: /cluster/work/climate/bestocke/data/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata 
#Data in my desktop is from: /Users/yunpeng/data/gpp_gmd/stocker20gmd_outputs/rdata_objects/out_eval_FULL.Rdata 
obs <-out_eval_FULL$gpp$fluxnet$data$meandf$obs
pred <- out_eval_FULL$gpp$fluxnet$data$meandf$mod
obs_pred <- as.data.frame(cbind(obs,pred))
obs_pred <- na.omit(obs_pred)
summary(obs_pred)

# According to textbook:Uncertainty estimates obtained as standard deviations of repeated measurement results are called A type uncertainty estimates. 
# In this way,  “sample mean” corresponds to the true GPP, and the “observation” is the modelled GPP.
obs_pred$variance <- (obs_pred$obs - obs_pred$pred)^2
uncertainty_gpp <- sqrt(sum(obs_pred$variance)/nrow(obs_pred))
uncertainty_gpp
#this is uncertainty of gpp
mean_gpp <- mean(obs_pred$pred,na.rm=TRUE)
gpp_pft <- gpp_df$gpp

# Uncertainty of TNPP/GPP
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

#now, calculate N uptake in the leaf = GPP * (ANPP/GPP) * (leafNPP/ANPP) * (leaf n/c) * (1-NRE)
uncertainty_lnf <- lnf_f * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                   (mod_anpp_uncertainty/(anpp_f/gpp_pft))^2 +
                                   (mod_lnpp_uncertainty/(lnpp_f/anpp_f))^2 +
                                   (mod_leafnc_uncertainty/leafnc_f)^2 +
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
summary(uncertainty_bnpp/bnpp_f)

uncertainty_bnf <- (bnf_f/bnpp_f) * uncertainty_bnpp 
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

#npp = 0.435 * gpp
uncertainty_grass_npp <- npp_g * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                         (summary(mod_tnpp_grass)$coef[1,2]/summary(mod_tnpp_grass)$coef[1,1])^2)

uncertainty_grass_npp_gpp <- summary(mod_tnpp_grass)$coef[1,2]

uncertainty_grass_anpp_gpp <- summary(mod_anpp_grass)$coef[1,2]



#anpp = 0.228 * gpp
uncertainty_grass_anpp <- anpp_g * sqrt( (uncertainty_gpp/mean_gpp)^2 +
                                           (summary(mod_anpp_grass)$coef[1,2]/summary(mod_anpp_grass)$coef[1,1])^2)

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


#now, output uncertainty values in table s1
#uncertainty for forest
sum(uncertainty_npp*(forest_percent *conversion)*available_grid2,na.rm=TRUE) #npp
sum(uncertainty_anpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE) #anpp
sum(uncertainty_bnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE) #bnpp
sum(uncertainty_lnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE) #leaf npp
sum(mod_lnpp_uncertainty*(forest_percent *conversion)*available_grid2,na.rm=TRUE)#leaf npp/gpp
sum(uncertainty_wnpp_gpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)#wood npp/gpp
sum(uncertainty_wnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)#wood npp
sum(uncertainty_lnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)#leaf n flux
sum(uncertainty_bnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)#root n flux
sum(uncertainty_wnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)#wood n flux
sum(uncertainty_nuptake*(forest_percent *conversion)*available_grid2,na.rm=TRUE)#n uptake 

#grassland
sum(uncertainty_grass_npp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)# npp
sum(uncertainty_grass_anpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)#anpp
sum(uncertainty_grass_bnpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)#bnpp
sum(uncertainty_grass_lnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)# leaf n flux
sum(uncertainty_grass_bnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)#root n flux
sum(uncertainty_grass_nuptake*(grass_percent *conversion)*available_grid2,na.rm=TRUE)# n uptake

#uncertainty of two pfts (forest + grassland)
# gpp
uncertainty_gpp
# npp/gpp
sqrt(sum(mod_tnpp_uncertainty*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + 
       sum(uncertainty_grass_npp_gpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
# npp
sqrt(sum(uncertainty_npp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_npp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
# anpp/gpp
sqrt(sum(mod_anpp_uncertainty*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + 
       sum(uncertainty_grass_anpp_gpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
# anpp
sqrt(sum(uncertainty_anpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_anpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
# bnpp
sqrt(sum(uncertainty_bnpp*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_bnpp*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
# leaf n flux
sqrt(sum(uncertainty_lnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_lnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
# root n flux
sqrt(sum(uncertainty_bnf*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_bnf*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
# n uptake
sqrt(sum(uncertainty_nuptake*(forest_percent *conversion)*available_grid2,na.rm=TRUE)^2 + sum(uncertainty_grass_nuptake*(grass_percent *conversion)*available_grid2,na.rm=TRUE)^2)
