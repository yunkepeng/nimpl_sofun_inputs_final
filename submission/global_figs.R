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
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")
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
load(file = "~/data/NPP_Grassland_final/statistical_model/tnpp_grass.RData")
mod_tnpp_grass<- tnpp_grass
summary(mod_tnpp_grass)
load(file = "~/data/NPP_Grassland_final/statistical_model/anpp_grass.RData")
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

#now, produce results
#now, produce results --> not used  plot_map4 at the end since its width is too larger shwon in map - not looks nice - if used then see comments below for example:
#devtools::load_all("/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/vhs/")
#devtools::load_all("/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/rbeni/")
# plot_map4(na.omit(all_maps[,c("lon","lat","gpp")]),varnam = "gpp")

#Figure 3 - map output along with validated sites in forest (red) and grassland (blue)
#1. gpp
NPP_forest <- read.csv("~/data/NPP_final/NPP_validation.csv")
NPP_grassland <- read.csv("~/data/NPP_Grassland_final/NPP_grass_validation.csv")

gpp_f <- (NPP_forest %>% filter(GPP>0) %>% filter(pred_gpp_c3>0))[,c("lon","lat")]
gpp_g <- (NPP_grassland %>% filter(GPP>0) %>% filter(weightedgpp_measured_c3>0))[,c("lon","lat")]
total_value <- round(sum(all_maps[,"gpp"]*conversion,na.rm=TRUE),2)

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","gpp")]),
                varnam = "gpp",latmin = -65, latmax = 85,combine=FALSE)
a1 <- gg$ggmap +
  geom_point(data=gpp_f,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=gpp_g,aes(lon,lat),col="blue",size=1.5)+
  labs(title = paste("GPP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)

a2 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))

#2. npp
npp_f <- (NPP_forest %>% filter(TNPP_1>0) %>% filter(pred_npp>0))[,c("lon","lat")]
npp_g <- (NPP_grassland %>% filter(TNPP_1>0) %>% filter(pred_npp>0))[,c("lon","lat")]
total_value <- round(sum(all_maps[,"npp_pft"]*conversion,na.rm=TRUE),2)

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","npp_pft")]),
                varnam = "npp_pft",latmin = -65, latmax = 85,combine=FALSE)

a3 <- gg$ggmap +
  geom_point(data=npp_f,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=npp_g,aes(lon,lat),col="blue",size=1.5)+
  labs(title = paste("BP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)

a4 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))

#3.anpp
anpp_f <- (NPP_forest %>% filter(ANPP_2>0) %>% filter(pred_anpp>0))[,c("lon","lat")]
anpp_g <- (NPP_grassland %>% filter(ANPP_2>0) %>% filter(pred_anpp>0))[,c("lon","lat")]
total_value <- round(sum(all_maps[,"anpp_pft"]*conversion,na.rm=TRUE),2)

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","anpp_pft")]),
                varnam = "anpp_pft",latmin = -65, latmax = 85,combine=FALSE)

a5 <- gg$ggmap +
  geom_point(data=anpp_f,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=anpp_g,aes(lon,lat),col="blue",size=1.5)+
  labs(title = paste("ANPP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)

a6 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))

#4. leaf c/n
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one 
SP_input2 <- SP_input[,c("lat","lon","z","Vcmax25","narea","lma")]
sitemean <- aggregate(SP_input2,by=list(SP_input2$lon,SP_input2$lat), FUN=mean, na.rm=TRUE) 
sitemean$pred_leafn <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax25/sitemean$lma
sitemean$obs_leafn <- sitemean$narea/sitemean$lma

laefcn <- (sitemean %>% filter(pred_leafn>0) %>% filter(obs_leafn>0))[,c("lon","lat")]

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","leafcn_pft")]),
                varnam = "leafcn_pft",latmin = -65, latmax = 85,combine=FALSE)

total_value <- round(mean(leafcn_pft,na.rm=TRUE),2)
a7 <- gg$ggmap +
  geom_point(data=laefcn,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("Leaf C/N: ", total_value))+
  theme_grey(base_size = 12)

a8 <- gg$gglegend

#5. NRE
NRE_validation <- read.csv(file="~/data/NPP_final/NRE_validation.csv") 
nre_site <- (NRE_validation %>% filter(pred_nre>0) %>% filter(NRE>0))[,c("lon","lat")]

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nre_pft")]),
                varnam = "nre_pft",latmin = -65, latmax = 85,combine=FALSE)

total_value <- round(mean(nre_pft,na.rm=TRUE),2)

a9 <- gg$ggmap +
  geom_point(data=nre_site,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("NRE: ", total_value))+
  theme_grey(base_size = 12)

a10 <- gg$gglegend

#6. nuptake
Nmin_validation <- read.csv("~/data/NPP_final/Nmin_validation.csv")
nuptake_f <- (Nmin_validation %>% filter(pred_nuptake>0) %>% filter(obs_nuptake>0))[,c("lon","lat")]

total_value <- 1000*round(sum(all_maps[,"nuptake_pft"]*conversion,na.rm=TRUE),2) #unit convert from PgN/yr to TgN/yr

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nuptake_pft")]),
                varnam = "nuptake_pft",latmin = -65, latmax = 85,combine=FALSE)

a11 <- gg$ggmap +
  geom_point(data=nuptake_f,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("N uptake: ", total_value, "TgN/yr", sep=" " ))+
  theme_grey(base_size = 12)

a12 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))



plot_grid(a1,a2,a3,a4,a5,a6,
          a7,a8,a9,a10,a11,a12,nrow=2,rel_widths = c(3/16, 1/16,3/16,1/16,3/16,1/16),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' '))

ggsave(paste("~/data/output/fig3.jpg",sep=""),width = 20, height = 10)

###Figure 5: CUE, NUE, Nuptake and its relation
all_maps$CUE <- all_maps$npp_pft/all_maps$gpp
total_value_CUE <- round(mean(all_maps$CUE,na.rm=TRUE),2)

all_maps$NUE <- all_maps$npp_pft/all_maps$nuptake_pft
total_value_NUE <- round(mean(all_maps$NUE,na.rm=TRUE),2)
#CUE
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","CUE")]),
                varnam = "CUE",latmin = -65, latmax = 85,combine=FALSE,
                breaks = seq(0.25,0.6,0.05))

b1 <- gg$ggmap +
  labs(title = paste("BPE: ", total_value_CUE))+
  theme_grey(base_size = 15)

b2 <- gg$gglegend
#NUE
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","NUE")]),
                varnam = "NUE",latmin = -65, latmax = 85,combine=FALSE,
                breaks = seq(45,85,5))

b3 <- gg$ggmap +
  labs(title = paste("NUE: ", total_value_NUE,"gC/gN"))+
  theme_grey(base_size = 15)

b4 <- gg$gglegend
#CUE~NUE
b5 <- ggplot(data=all_maps, aes(x=CUE, y=NUE)) +
  geom_point(aes(color=nuptake_pft),alpha=0.3,size=0.3)+geom_smooth(method = "lm", se = TRUE)+ 
  scale_color_viridis(discrete=FALSE,direction= -1)+theme_classic()+theme(axis.title = element_text(size = 30),
                                                                          axis.text = element_text(size = 20),
                                                                          legend.title = element_text(size = 14))+
  xlab("BPE")+ylab("NUE (gC/gN)")+labs(color= ~paste("N uptake", " (gN m"^-2,"yr"^-1,")"))
summary(lm(all_maps$NUE~all_maps$CUE))

plot_grid(b1,b2,b3,b4,b5,
          labels = c('(a)',' ','(b)',' ','(c)'), rel_widths = c(3/12, 1/21,3/12,1/12,4/12),
          nrow=1,label_size = 20)

ggsave(paste("~/data/output/fig5.jpg",sep=""),width = 20, height = 5)
#variation calculation
#calculate nup = gpp * bpe /nue, about variation of each component
a1 <- log((all_maps$gpp * all_maps$CUE * all_maps$NUE)/(mean(all_maps$gpp,na.rm = TRUE) * all_maps$CUE * all_maps$NUE)) #gpp
a2 <- log((all_maps$gpp * all_maps$CUE * all_maps$NUE)/(all_maps$gpp* mean(all_maps$CUE,na.rm = TRUE) * all_maps$NUE)) #cue
a3 <- log((all_maps$gpp * all_maps$CUE * all_maps$NUE)/(all_maps$gpp* all_maps$CUE * mean(all_maps$NUE,na.rm = TRUE))) #nue
a_all <- as.data.frame(cbind(a1,a2,a3))
total_sum <- sum(abs((a_all)*conversion),na.rm=TRUE)
sum(abs((a1)*conversion),na.rm=TRUE)/total_sum
sum(abs((a2)*conversion),na.rm=TRUE)/total_sum
sum(abs((a3)*conversion),na.rm=TRUE)/total_sum


#computating correltion coefficient r
cor(na.omit(all_maps$CUE), na.omit(all_maps$NUE),
    method = c("pearson", "kendall", "spearman"))

#Fig. S3 N uptake and NUE in forest and grassland
#Nuptake
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nuptake_pft")]),
                varnam = "nuptake_pft",latmin = -65, latmax = 85,combine=FALSE)

c1 <- gg$ggmap +labs(title = "N uptake")+theme_grey(base_size = 20)

c2 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nuptake_forest")]),
                varnam = "nuptake_forest",latmin = -65, latmax = 85,combine=FALSE)

c3 <- gg$ggmap +labs(title = "Forest N uptake")+theme_grey(base_size = 20)

c4 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))

c5 <- ggplot(data=all_maps, aes(all_maps$nuptake_pft)) + 
  geom_histogram()+
  geom_vline(aes(xintercept = mean(all_maps$nuptake_forest,na.rm=TRUE),col='forest'),size=2)+
  geom_vline(aes(xintercept = mean(all_maps$nuptake_grass,na.rm=TRUE),col='grassland'),size=2)+
  theme(text = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))+
  theme(legend.title = element_blank())+
  labs(x = ~paste("Total N uptake (","gN m"^-2,"yr"^-1, ")"))

#NUE
all_maps$NUE_forest <- all_maps$npp_forest/all_maps$nuptake_forest
all_maps$NUE_grass <- all_maps$npp_grass/all_maps$nuptake_grass

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","NUE")]),
                varnam = "NUE",latmin = -65, latmax = 85,combine=FALSE,
                breaks = seq(45,90,5))

c6 <- gg$ggmap +labs(title = "NUE (gC/gN)")+theme_grey(base_size = 20)

c7 <- gg$gglegend

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","NUE_forest")]),
                varnam = "NUE_forest",latmin = -65, latmax = 85,combine=FALSE,
                breaks = seq(45,90,5))

c8 <- gg$ggmap +labs(title = "Forest NUE (gC/gN)")+theme_grey(base_size = 20)

c9 <- gg$gglegend

c10 <- ggplot(data=all_maps, aes(all_maps$NUE)) + 
  geom_histogram()+
  geom_vline(aes(xintercept = mean(all_maps$NUE_forest,na.rm=TRUE),col='forest'),size=2)+
  geom_vline(aes(xintercept = mean(all_maps$NUE_grass,na.rm=TRUE),col='grassland'),size=2)+
  theme(text = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))+
  theme(legend.title = element_blank())+
  labs(x = ~paste("NUE (gC/gN)"))

plot_grid(c1,c2,c3,c4,c5,
          c6,c7,c8,c9,c10,nrow=2,rel_widths = c(3/11, 1/11,3/11,1/11,3/11),
          labels = c('(a)',' ','(b)',' ','(c)',
                     '(d)',' ','(e)',' ','(f)'),label_size = 20)

ggsave(paste("~/data/output/figS3.jpg",sep=""),width = 20, height = 10)

#fig.s4 - representing all predictors
gg <- plot_map3(na.omit(CNrt[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15);d2 <- gg$gglegend

gg <- plot_map3(na.omit(age[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15);d4 <- gg$gglegend+labs(title = ~paste("years"))

gg <- plot_map3(na.omit(fAPAR[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15);d6 <- gg$gglegend

gg <- plot_map3(na.omit(PPFD[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE,breaks = seq(100,800,100))
d7 <- gg$ggmap +labs(title = "PPFD")+theme_grey(base_size = 15);d8 <- gg$gglegend + labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(Tg[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d9 <- gg$ggmap +labs(title =~paste(T[g]))+theme_grey(base_size = 15);d10 <- gg$gglegend + labs(title =~paste("\u00B0C"))

gg <- plot_map3(na.omit(vpd[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d11 <- gg$ggmap +labs(title = "D")+theme_grey(base_size = 15);d12 <- gg$gglegend+ labs(title =~paste("kPa"))

gg <- plot_map3(na.omit(vcmax25_df[,c("lon","lat","vcmax25")]), varnam = "vcmax25",latmin = -65, latmax = 85,combine=FALSE)
d13 <- gg$ggmap +labs(title =~paste(V[cmax25]))+theme_grey(base_size = 15);d14 <- gg$gglegend+ labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(LMA[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d15 <- gg$ggmap +labs(title = "LMA")+theme_grey(base_size = 15);d16 <- gg$gglegend+ labs(title =~paste("g ","m"^-2))

plot_grid(d1,d2,d3,d4,d5,d6,
          d7,d8,d9,d10,d11,d12,
          d13,d14,d15,d16,
          nrow=3,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' ',
                     '(g)',' ','(h)',' '),label_size = 15)

ggsave(paste("~/data/output/figS4.jpg",sep=""),width = 20, height = 10)

#now - moving to final part - environmental factors on CUE/NUE (Fig.4), CUE(Fig.S1) and NUE(Fig.S2)

#create a function to control each factor --> output ratio of N demand of GPP
nuptake_pft_final <- all_maps$nuptake_pft # this map is by 'default' --> since our calculation accounts for log( Nup/GPP / Nup`/GPP), this is Nup
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
  
  leafnc_f <- (summary(n1)$coef[1,1]/0.46) + 
    (summary(n1)$coef[2,1]/0.46) *vcmax25_pred/LMA_pred
  #0.46 is constant Cmass
  
  nre_f <- (1/(1+exp(-(summary(nre_model)$coef[1,1]+
                         summary(nre_model)$coef[2,1] *Tg_pred + 
                         summary(nre_model)$coef[3,1] * log(vpd_pred)))))
  
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
  
  #pft
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
  
  nre_pft <- available_grid2*(nre_f*forest_percent +nre_g*grass_percent)
  nre_forest <-  available_grid2*nre_f*forest_percent
  nre_grassland <- available_grid2*nre_g*grass_percent
  
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
  nuptake_pft_ratio <- log((nuptake_pft_final/gpp_df$gpp)/(nuptake_pft/gpp_df$gpp))
  return(nuptake_pft_ratio)
}

nuptake_standard <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
summary(nuptake_standard) # now grid is 0 since nothing changed now --> Nup = Nup` there log(1) = 0 for all grids

#firstly, get mean value of each predictors
length(Tg$myvar[is.na(available_grid2)==FALSE]) # 55933 grids are possible for all predictors
mean_Tg <- mean(Tg$myvar[is.na(available_grid2)==FALSE]);mean_Tg
mean_PPFD <- mean(PPFD$myvar[is.na(available_grid2)==FALSE]);mean_PPFD
mean_vpd <- mean(vpd$myvar[is.na(available_grid2)==FALSE]);mean_vpd
mean_fAPAR <- mean(fAPAR$myvar[is.na(available_grid2)==FALSE]);mean_fAPAR
mean_age <- mean(age$myvar[is.na(available_grid2)==FALSE]);mean_age
mean_CNrt <- mean(CNrt$myvar[is.na(available_grid2)==FALSE]);mean_CNrt
mean_LMA <- mean(LMA$myvar[is.na(available_grid2)==FALSE]);mean_LMA
mean_vcmax25 <- mean(vcmax25_df$vcmax25[is.na(available_grid2)==FALSE]);mean_vcmax25

nuptake_Tg <- cal_nuptake(rep(mean_Tg,259200),PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_PPFD <- cal_nuptake(Tg$myvar,rep(mean_PPFD,259200),vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_vpd <- cal_nuptake(Tg$myvar,PPFD$myvar,rep(mean_vpd,259200),fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_fAPAR <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,rep(mean_fAPAR,259200),age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_age <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,rep(mean_age,259200),CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nuptake_CNrt <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,rep(mean_CNrt,259200),LMA$myvar,vcmax25_df$vcmax25)
nuptake_LMA <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,rep(mean_LMA,259200),vcmax25_df$vcmax25)
nuptake_vcmax25 <- cal_nuptake(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,rep(mean_vcmax25,259200))

nuptake_all <- as.data.frame(cbind(vcmax25_df,nuptake_Tg,nuptake_PPFD,nuptake_vpd,nuptake_fAPAR,nuptake_age,nuptake_CNrt,nuptake_LMA,nuptake_vcmax25))
summary(nuptake_all)

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_CNrt")]),varnam = "nuptake_CNrt",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.2,0.2,0.05))
e1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15);e2 <- gg$gglegend

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_age")]),varnam = "nuptake_age",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3","tomato4"),
                breaks = seq(-0.20,0.40,0.050))
e3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15);e4 <- gg$gglegend+labs(title = ~paste("years"))

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_fAPAR")]),varnam = "nuptake_fAPAR",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3","tomato4"),
                breaks = seq(-0.20,0.40,0.050))
e5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15);e6 <- gg$gglegend

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_PPFD")]),varnam = "nuptake_PPFD",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
e7 <- gg$ggmap +labs(title = "PPFD")+theme_grey(base_size = 15);e8 <- gg$gglegend + labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_Tg")]),varnam = "nuptake_Tg",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.40,0.40,0.10))
e9 <- gg$ggmap +labs(title =~paste(T[g]))+theme_grey(base_size = 15);e10 <- gg$gglegend + labs(title =~paste("\u00B0C"))

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_vpd")]),varnam = "nuptake_vpd",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.20,0.20,0.05))
e11 <- gg$ggmap +labs(title = "D")+theme_grey(base_size = 15);e12 <- gg$gglegend+ labs(title =~paste("kPa"))

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_vcmax25")]), varnam = "nuptake_vcmax25",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.20,0.20,0.05))
e13 <- gg$ggmap +labs(title =~paste(V[cmax25]))+theme_grey(base_size = 15);e14 <- gg$gglegend+ labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nuptake_all[,c("lon","lat","nuptake_LMA")]),varnam = "nuptake_LMA",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
e15 <- gg$ggmap +labs(title = "LMA")+theme_grey(base_size = 15);e16 <- gg$gglegend+ labs(title =~paste("g ","m"^-2))

plot_grid(e1,e2,e3,e4,e5,e6,
          e7,e8,e9,e10,e11,e12,
          e13,e14,e15,e16,
          nrow=3,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' ',
                     '(g)',' ','(h)',' '),label_size = 15)

ggsave(paste("~/data/output/fig4.jpg",sep=""),width = 20, height = 10)

#now, work on environmental factors on CUE (fig.s1)
cue_final <- all_maps$CUE # by default
cal_cue <- function(fAPAR_pred,age_pred,CNrt_pred){
  #forest
  npp_f <- gpp_df$gpp * (1/(1 + exp(-(summary(mod_tnpp)$coef[1,1]+
                                        summary(mod_tnpp)$coef[2,1]* log(CNrt_pred)+
                                        summary(mod_tnpp)$coef[3,1] * log(age_pred) + 
                                        summary(mod_tnpp)$coef[4,1]* fAPAR_pred))))
  #grass
  npp_g <- gpp_df$gpp * summary(mod_tnpp_grass)$coef[1,1]
  
  #pft
  npp_pft <- available_grid2* (npp_f*forest_percent +npp_g*grass_percent)
  cue_pft <- npp_pft/gpp_df$gpp
  cue_pft_ratio <- log(cue_final/cue_pft)
  return(cue_pft_ratio)
}

cue_standard <- cal_cue(fAPAR$myvar,age$myvar,CNrt$myvar)
summary(cue_standard) # by default
cue_fAPAR <- cal_cue(rep(mean_fAPAR,259200),age$myvar,CNrt$myvar)
cue_age <- cal_cue(fAPAR$myvar,rep(mean_age,259200),CNrt$myvar)
cue_CNrt <- cal_cue(fAPAR$myvar,age$myvar,rep(mean_CNrt,259200))
cue_all <- as.data.frame(cbind(vcmax25_df,cue_fAPAR,cue_age,cue_CNrt))

gg <- plot_map3(na.omit(cue_all[,c("lon","lat","cue_CNrt")]),varnam = "cue_CNrt",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.1,0.1,0.02))
f1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15);f2 <- gg$gglegend

gg <- plot_map3(na.omit(cue_all[,c("lon","lat","cue_age")]),varnam = "cue_age",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3","tomato4"),
                breaks = seq(-0.20,0.40,0.050))
f3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15);f4 <- gg$gglegend+labs(title = ~paste("years"))

gg <- plot_map3(na.omit(cue_all[,c("lon","lat","cue_fAPAR")]),varnam = "cue_fAPAR",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3","tomato4"),
                breaks = seq(-0.20,0.40,0.050))
f5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15);f6 <- gg$gglegend

plot_grid(f1,f2,f3,f4,f5,f6,
          nrow=1,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' '),label_size = 15)

ggsave(paste("~/data/output/figS1.jpg",sep=""),width = 20, height = 5)

#now, work on nue
nue_final <- all_maps$NUE # by default
cal_nue <- function(Tg_pred,PPFD_pred,vpd_pred,fAPAR_pred,age_pred,CNrt_pred,LMA_pred,vcmax25_pred){
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
  
  leafnc_f <- (summary(n1)$coef[1,1]/0.46) + 
    (summary(n1)$coef[2,1]/0.46) *vcmax25_pred/LMA_pred
  #0.46 is constant Cmass
  
  nre_f <- (1/(1+exp(-(summary(nre_model)$coef[1,1]+
                         summary(nre_model)$coef[2,1] *Tg_pred + 
                         summary(nre_model)$coef[3,1] * log(vpd_pred)))))
  
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
  
  #pft
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
  
  nre_pft <- available_grid2*(nre_f*forest_percent +nre_g*grass_percent)
  nre_forest <-  available_grid2*nre_f*forest_percent
  nre_grassland <- available_grid2*nre_g*grass_percent
  
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
  nue_pft_ratio <- log(nue_final/(npp_pft/nuptake_pft))
  return(nue_pft_ratio)
}

nue_standard <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
summary(nue_standard)

nue_Tg <- cal_nue(rep(mean_Tg,259200),PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nue_PPFD <- cal_nue(Tg$myvar,rep(mean_PPFD,259200),vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nue_vpd <- cal_nue(Tg$myvar,PPFD$myvar,rep(mean_vpd,259200),fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nue_fAPAR <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,rep(mean_fAPAR,259200),age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nue_age <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,rep(mean_age,259200),CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25)
nue_CNrt <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,rep(mean_CNrt,259200),LMA$myvar,vcmax25_df$vcmax25)
nue_LMA <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,rep(mean_LMA,259200),vcmax25_df$vcmax25)
nue_vcmax25 <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,rep(mean_vcmax25,259200))

nue_all <- as.data.frame(cbind(vcmax25_df,nue_Tg,nue_PPFD,nue_vpd,nue_fAPAR,nue_age,nue_CNrt,nue_LMA,nue_vcmax25))
summary(nue_all)

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_CNrt")]),varnam = "nue_CNrt",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.02,0.02,0.005))
g1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15);g2 <- gg$gglegend

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_age")]),varnam = "nue_age",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4","royalblue3","wheat","tomato3"),
                breaks = seq(-0.06,0.04,0.01))
g3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15);g4 <- gg$gglegend+labs(title = ~paste("years"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_fAPAR")]),varnam = "nue_fAPAR",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4","royalblue3", "wheat","tomato3"),
                breaks = seq(-0.06,0.04,0.01))
g5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15);g6 <- gg$gglegend

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_PPFD")]),varnam = "nue_PPFD",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3","tomato4"),
                breaks = seq(-0.08,0.12,0.02))
g7 <- gg$ggmap +labs(title = "PPFD")+theme_grey(base_size = 15);g8 <- gg$gglegend + labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_Tg")]),varnam = "nue_Tg",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.40,0.40,0.10))
g9 <- gg$ggmap +labs(title =~paste(T[g]))+theme_grey(base_size = 15);g10 <- gg$gglegend + labs(title =~paste("\u00B0C"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_vpd")]),varnam = "nue_vpd",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.20,0.20,0.05))
g11 <- gg$ggmap +labs(title = "D")+theme_grey(base_size = 15);g12 <- gg$gglegend+ labs(title =~paste("kPa"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_vcmax25")]), varnam = "nue_vcmax25",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.20,0.20,0.05))
g13 <- gg$ggmap +labs(title =~paste(V[cmax25]))+theme_grey(base_size = 15);g14 <- gg$gglegend+ labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_LMA")]),varnam = "nue_LMA",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
g15 <- gg$ggmap +labs(title = "LMA")+theme_grey(base_size = 15);g16 <- gg$gglegend+ labs(title =~paste("g ","m"^-2))

plot_grid(g1,g2,g3,g4,g5,g6,
          g7,g8,g9,g10,g11,g12,
          g13,g14,g15,g16,
          nrow=3,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' ',
                     '(g)',' ','(h)',' '),label_size = 15)

ggsave(paste("~/data/output/figS2.jpg",sep=""),width = 20, height = 10)

#table S1 - global estimations and uncertainty evaluation
#####Table S1
#just print estimations of all values
for (i in 4:31){
  varname <- names(all_maps)[i]
  total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
  print(varname)
  print(total_value)
}