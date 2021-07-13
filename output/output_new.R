rm(list=ls())
#for figure 1 - just summarizing each part

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

leafcn_pft <- 1/(available_grid2*leafcn_forest * forest_percent+
  available_grid2* (1/18) * grass_percent)

nre_pft <- available_grid2*nre_f * forest_percent+
  available_grid2*0.69 * grass_percent

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
                                lnpp_forest,wnpp_forest,leafcn_forest,nre_pft,
                                leafcn_pft,wnf_forest,
                                lnf_pft,lnf_forest,lnf_grass,
                                bnf_pft,bnf_forest,bnf_grass,
                                nuptake_pft,nuptake_forest,nuptake_grass))
all_maps$cue_pft <- all_maps$npp_pft/all_maps$gpp
all_maps$nue_pft <- all_maps$npp_pft/all_maps$nuptake_pft

plot_map3(all_maps[,c("lon","lat","npp_pft")],
          varnam = "npp_pft",
          latmin = -65, latmax = 85,font_size = 12)
ggsave(paste("/Users/yunpeng/data/output/output_map/npp_pft.jpg",sep=""))

plot_map3(all_maps[,c("lon","lat","leafcn_pft")],
          varnam = "leafcn_pft",
          latmin = -65, latmax = 85,font_size = 12)
ggsave(paste("/Users/yunpeng/data/output/output_map/leafcn_pft.jpg",sep=""))

plot_map3(all_maps[,c("lon","lat","nre_pft")],
          varnam = "nre_pft",
          latmin = -65, latmax = 85,font_size = 12)
ggsave(paste("/Users/yunpeng/data/output/output_map/nre_pft.jpg",sep=""))

plot_map3(all_maps[,c("lon","lat","cue_pft")],
          varnam = "cue_pft",
          latmin = -65, latmax = 85,font_size = 12)
ggsave(paste("/Users/yunpeng/data/output/output_map/cue_pft.jpg",sep=""))

plot_map3(all_maps[,c("lon","lat","nue_pft")],
          varnam = "nue_pft",
          latmin = -65, latmax = 85,font_size = 12)
ggsave(paste("/Users/yunpeng/data/output/output_map/nue_pft.jpg",sep=""))

plot_map3(all_maps[,c("lon","lat","nuptake_pft")],
          varnam = "nuptake_pft",
          latmin = -65, latmax = 85,font_size = 12)
ggsave(paste("/Users/yunpeng/data/output/output_map/nuptake_pft.jpg",sep=""))


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
GPP <- gpp_df$gpp
CUE <- npp_pft/gpp_df$gpp
NUE <- npp_pft/nuptake_pft
Nup <- GPP*CUE*(1/NUE)
#test
summary(log(Nup/(GPP*CUE*(1/NUE))))

n_GPP <- log(Nup/(mean(GPP,na.rm=TRUE)*CUE*(1/NUE)))
n_CUE <- log(Nup/(GPP*mean(CUE,na.rm=TRUE)*(1/NUE)))
n_NUE <- log(Nup/(GPP*CUE*(1/mean(NUE,na.rm=TRUE))))
n_all <- as.data.frame(cbind(n_GPP,n_CUE,n_NUE))
total_sum <- sum(abs((n_all)*conversion),na.rm=TRUE)
#now, research on 4 componenets
sum(abs((n_GPP)*conversion),na.rm=TRUE)/total_sum
sum(abs((n_CUE)*conversion),na.rm=TRUE)/total_sum
sum(abs((n_NUE)*conversion),na.rm=TRUE)/total_sum

all_maps1 <- as.data.frame(cbind(gpp_df,n_GPP,n_CUE,n_NUE))
summary(all_maps1)
plot_map3(all_maps1[,c("lon","lat","n_GPP")],
          varnam = "n_GPP",
          latmin = -65, latmax = 85,font_size = 15,
          colorscale = c( "royalblue4", "wheat","tomato3"),
          breaks = seq(-3,3,0.5))

ggsave(paste("/Users/yunpeng/data/output/output_map/effect_GPP.jpg",sep=""))

plot_map3(all_maps1[,c("lon","lat","n_CUE")],
          varnam = "n_CUE",
          latmin = -65, latmax = 85,font_size = 15,
          colorscale = c( "royalblue4", "wheat","tomato3"),
          breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1, 0,
                     0.1,0.2,0.3,0.4,0.5))
ggsave(paste("/Users/yunpeng/data/output/output_map/effect_n_CUE.jpg",sep=""))

plot_map3(all_maps1[,c("lon","lat","n_NUE")],
          varnam = "n_NUE",
          latmin = -65, latmax = 85,font_size = 15,
          colorscale = c( "royalblue4", "wheat","tomato3"),
          breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1, 0,
                     0.1,0.2,0.3,0.4,0.5))
ggsave(paste("/Users/yunpeng/data/output/output_map/effect_n_NUE.jpg",sep=""))

summary(all_maps)
My_Theme = theme(
  axis.title.x = element_text(size = 40),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 40),
  axis.text.y = element_text(size = 20))

ggplot(data=all_maps, aes(x=cue_pft, y=nue_pft)) +
  geom_point()+geom_smooth(method = "lm", se = TRUE)+
  xlab("CUE")+ylab("NUE")+theme_classic()+My_Theme

mid <- mean(all_maps$nuptake_pft,na.rm=TRUE)

ggplot(data=all_maps, aes(x=cue_pft, y=nue_pft)) +
  geom_point(aes(color=nuptake_pft),alpha=0.5)+geom_smooth(method = "lm", se = TRUE)+ 
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "white",high = "red", space = "Lab")+
  xlab("CUE")+ylab("NUE")+theme_classic()+theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 40),
    axis.text.y = element_text(size = 20))

ggsave(paste("/Users/yunpeng/data/output/output_map/NUE_CUE.jpg",sep=""))
summary(lm(nue_pft~cue_pft,all_maps))

all_maps1_max <- all_maps1[,5:8]

for (i in 1:nrow(all_maps1_max)){
  if (is.na(all_maps1_max[i,1])==TRUE){
    all_maps1_max$most_factor[i] <- NA
    all_maps1_max$most_factor_value[i] <- NA
  } else {
    all_maps1_max$most_factor[i] <- names((which.max(abs((all_maps1_max[i,1:4])))))
    all_maps1_max$most_factor_value[i] <- max(abs((all_maps1_max[i,1:4])))
  }
}

apparent_point <- as.data.frame(cbind(gpp_df[,c("lon","lat")],all_maps1_max[,c("most_factor","most_factor_value")]))
apparent_point_available <- subset(apparent_point,is.na(most_factor_value)==FALSE) 
dim(apparent_point_available)
apparent_point_available$most_factor_value <- as.numeric(apparent_point_available$most_factor_value)

apparent_point_available %>% group_by(most_factor) %>% summarise(number=n())

area_final <- as.data.frame(cbind(gpp_df[,c("lon","lat")],all_maps1_max[,c("most_factor","most_factor_value")]))
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
                varnam = "nuptake_pft",plot_title = paste(" "),
                latmin = -65, latmax = 85,combine=FALSE)
apparent_point_available$most_factor[apparent_point_available$most_factor=="n_CUE"] <- "CUE"
apparent_point_available$most_factor[apparent_point_available$most_factor=="n_LUE"] <- "LUE"
apparent_point_available$most_factor[apparent_point_available$most_factor=="n_NUE_1"] <- "N demand of NPP"
apparent_point_available$most_factor[apparent_point_available$most_factor=="n_PAR"] <- "PAR"

colors <-  c("red","red","red","red","red","red","red","red","red","black","yellow","green","orange","red","blue","purple")
gg$ggmap + geom_point(data=apparent_point_available,aes(lon,lat,color=most_factor),size=0.5)+
  scale_color_manual(values = colors)+ theme(
    legend.text = element_text(size = 20))+
  guides(colour = guide_legend(override.aes = list(size = 5)))
ggsave(paste("/Users/yunpeng/data/output/output_onefactor/allfactor_fig1.jpg",sep=""))



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
