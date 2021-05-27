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

save.image(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/uncertainty/uncertainty.Rdata")

load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/uncertainty/uncertainty.Rdata")

#######input GMD gpp to check its consistency with gmd gpp
#Data in Euler is from: /cluster/work/climate/bestocke/data/sofun_outputs/output_nc_global/global_FULL_fAPAR3g_v2_2000_2016.a.gpp.nc
#Data in my desktop is from: /Users/yunpeng/data/gpp_gmd/global_FULL_fAPAR3g_v2_2000_2016.a.gpp.nc
ncin <- nc_open(paste ("/Users/yunpeng/data/gpp_gmd/global_FULL_fAPAR3g_v2_2000_2016.a.gpp.nc"))
lon <- ncvar_get(ncin,"lon")
lat<-ncvar_get(ncin,"lat")
gpp <- ncvar_get(ncin,"gpp")
dim(gpp)
nc_close(ncin)
pre.vec.long <- as.vector(gpp)
pre.mat <- matrix(pre.vec.long, nrow = 259200, ncol = 17)
lonlat <- expand.grid(lon, lat)
gpp_gmd <- as.data.frame(cbind(lonlat,rowMeans(pre.mat)))
gpp_gmd$nimpl_gpp <- gpp_df$gpp
#compare nimpl gpp vs. gmd gpp
names(gpp_gmd) <- c("lon","lat","GMD_gpp","nimpl_gpp")
plot(GMD_gpp~nimpl_gpp,gpp_gmd)
analyse_modobs2(subset(gpp_gmd,nimpl_gpp>0),"GMD_gpp","nimpl_gpp", type = "points")

####check everything here
###forest
#check all output here...between fortran and R it should be all consistent (except for some values = 0 (NA) due to edge of grid in Fortran)
npp_df$npp_r <- gpp_df$gpp * (1/(1 + exp(-(-0.36075 * log(CNrt$myvar) -0.16213 * log(age$myvar) + 
                                             0.72793 * fAPAR$myvar + 0.57014))))
plot(npp_r~npp,subset(npp_df,npp>0))
anpp_df$anpp_r <- gpp_df$gpp * (1/(1 + exp(-(-0.55151 * log(CNrt$myvar) -0.20050 * log(age$myvar) + 
                                               1.06611 * fAPAR$myvar+ 0.35817))))

plot(anpp_r~anpp,subset(anpp_df,anpp>0))
aaa <- (subset(anpp_df,anpp==0&anpp_r!=0))
dim(aaa)
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(aaa$lon,aaa$lat, col="red", pch=16,cex=1)

bnpp_df$bnpp_r <- npp_df$npp_r- anpp_df$anpp_r
plot(bnpp_r~bnpp,subset(bnpp_df,bnpp>0))

lnpp_df$lnpp_r <- anpp_df$anpp_r* (1/(1 + exp(-(0.97093* log(PPFD$myvar) +
                                                  0.06453 * (Tg$myvar) 
                                                -0.80397 * log(vpd$myvar)
                                                -7.47165))))
plot(lnpp_r~lnpp,lnpp_df)
aaa <- (subset(lnpp_df,lnpp==0&lnpp_r!=0))
dim(aaa)
library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(aaa$lon,aaa$lat, col="red", pch=16,cex=1)


wnpp_df$wnpp_r <- anpp_df$anpp_r-lnpp_df$lnpp_r
plot(wnpp_r~wnpp,subset(wnpp_df,wnpp>0))

leafcn_df$leafcn_r <- (0.01599/0.46) + (0.005992/0.46) *vcmax25_df$vcmax25/LMA$myvar
plot(leafcn_r~leafcn,subset(leafcn_df,leafcn>0))

nre_df$nre_r <- (1/(1+exp(-(-0.064460 *Tg$myvar + 0.402850 * log(vpd$myvar) + 1.368935))))
plot(nre_r~nre,subset(nre_df,nre>0))

lnf_df$lnf_r <- (1-nre_df$nre_r)* leafcn_df$leafcn_r* lnpp_df$lnpp_r
plot(lnf_r~lnf,subset(lnf_df,lnf>0))

wnf_df$wnf_r <- wnpp_df$wnpp_r/100
plot(wnf_r~wnf,subset(wnf_df,wnf>0))

bnf_df$bnf_r <- bnpp_df$bnpp_r/94
plot(bnf_r~bnf,subset(bnf_df,bnf>0))

nuptake_df$nuptake_r <- lnf_df$lnf_r +wnf_df$wnf_r+bnf_df$bnf_r 
plot(nuptake_r~nuptake,subset(nuptake_df,nuptake>0))


###grassland
npp_grass_df$npp_r <- gpp_df$gpp *0.435
anpp_grass_df$anpp_r <- gpp_df$gpp *0.228
bnpp_grass_df$bnpp_r <- npp_grass_df$npp_r - anpp_grass_df$anpp_r 
lnf_grass_df$lnf_r <- anpp_grass_df$anpp_r *(1/18)*(1-(1/(1+exp(-(-0.064460 *Tg$myvar + 0.402850 * log(vpd$myvar) + 1.368935)))))
bnf_grass_df$bnf_r <- bnpp_grass_df$bnpp_r *(1/41)
nuptake_grass_df$nuptake_r <- lnf_grass_df$lnf_r + bnf_grass_df$bnf_r

plot(npp_grass~npp_r,npp_grass_df)
plot(anpp_grass~anpp_r,anpp_grass_df)
plot(bnpp_grass~bnpp_r,bnpp_grass_df)
plot(lnf_grass~lnf_r,subset(lnf_grass_df,lnf_grass>0))
plot(bnf_grass~bnf_r,bnf_grass_df)
plot(nuptake_grass~nuptake_r,subset(nuptake_grass_df,nuptake_grass>0))
summary(gpp_df)

# ALL consistent now!