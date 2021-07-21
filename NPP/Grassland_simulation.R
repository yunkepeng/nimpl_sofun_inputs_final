#re-processing fapar; input sites firstly
rm(list=ls())
library(dplyr)
library(tidyverse)  # depends
library(ncmeta)
library(viridis)
library(ggthemes)
library(LSD)
library(yardstick)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(tidyselect)
library(extrafont)
library(rbeni)
library(raster)
library(spgwr)
library(maps)
library(rworldmap)
library(cowplot)
library(spgwr)
library(lubridate)
library(lme4)
library(MuMIn)
library(lmerTest)
#load(file = "/Users/yunpeng/data/NPP_Grassland_final/grass_simulation.Rdata")
devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/rsofun/")
#see "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_site_orig.R" - sitename and sitename_fapar already been well checked
NPP_grassland  <- read.csv("/Users/yunpeng/data/NPP_Grassland_final/NPP_grassland.csv")

NPP_grassland$sitename <- NA
#not used this sitename (NPP1...) anymore - all
#climate and fapar both using sitemame2 (ingestr1-485)
for (i in 1:nrow(NPP_grassland)){
  NPP_grassland$sitename[i] <- paste("NPP",i,sep = "")
}


siteinfo <- data.frame(
  sitename = NPP_grassland$sitename,
  lon = NPP_grassland$lon,
  lat = NPP_grassland$lat,
  elv = NPP_grassland$z,
  year_start = NPP_grassland$Begin_year,
  year_end = NPP_grassland$End_year
)

siteinfo$no <- c(1:nrow(siteinfo))

siteinfo$year_start[siteinfo$year_start<=1980] <- 1980
siteinfo$year_end[siteinfo$year_start<=1980] <- 1989

siteinfo <-  siteinfo %>% dplyr::mutate(date_start = lubridate::ymd(paste0(year_start, "-01-01"))) %>%
  dplyr::mutate(date_end = lubridate::ymd(paste0(year_end, "-12-31"))) 


#aggregate by lon, lat, z, year_start, year_end
siteinfo2 <- aggregate(siteinfo,by=list(siteinfo$lon,siteinfo$lat,siteinfo$elv,siteinfo$year_start,siteinfo$year_end), FUN=mean, na.rm=TRUE) #site-mean
for (i in 1:nrow(siteinfo2)){
  siteinfo2$sitename2[i] <- paste("ingestr",i,sep = "")
}

siteinfo3 <- siteinfo2[,c("sitename2","lon","lat","elv","year_start","year_end")]

NPP_grassland_all <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat","elv","year_start","year_end"),all.x=TRUE), 
                           list(siteinfo,siteinfo3))

NPP_grassland_all <- NPP_grassland_all[order(NPP_grassland_all$no), ]
#we will use above data at the end
summary(NPP_grassland_all)

siteinfo_final <- siteinfo2[,c("sitename2","lon","lat","elv","year_start","year_end","date_start","date_end")]
names(siteinfo_final) <- c("sitename","lon","lat","elv","year_start","year_end","date_start","date_end")
summary(siteinfo_final)
dim(siteinfo_final) #485 sites to be obtained

#last check - if year_end and year all behave well
stopifnot( all(siteinfo_final$year_start == floor(siteinfo_final$year_start)) )
stopifnot( all(siteinfo_final$year_end == floor(siteinfo_final$year_end)) )
summary(siteinfo_final$year_end  - siteinfo_final$year_start)

#csvfile <- paste("/Users/yunpeng/data/NPP_final/fpar_name/grassland_fpar_name.csv",sep = "")
#write.csv(siteinfo_final, csvfile, row.names = TRUE)

#1. input fapar
fapar_df <- list.files("/Users/yunpeng/data/NPP_Grassland_final/reprocessing_fapar_final/fpar_all/",full.names = T)
length(fapar_df) # expected that 10 points were missing - due to no fapar data available in n_focal = 0, or 1, or 2. See reprocessing_fapar_all15yrs

for (i in 1:length(fapar_df)){
  df1 <- read.csv(fapar_df[i])
  df1$date <- as.Date(df1$date)
  df1$Year <- as.numeric(format(df1$date, format="%Y"))
  df1 <- df1[!(format(df1$date,"%m") == "02" & format(df1$date, "%d") == "29"), , drop = FALSE]
  site_name <- substr(sub('.*daily_', '', fapar_df[i]),1,nchar(sub('.*daily_', '', fapar_df[i]))-4) #specify site_name, measurement year of start and end
  yr_start <- subset(siteinfo_final,siteinfo_final$sitename == site_name)$year_start
  yr_end <- subset(siteinfo_final,siteinfo_final$sitename == site_name)$year_end
  
  if (yr_start<=2002) { # if measurement year before 2002 for a certain site --> calculating 2003-2012 average of MCD15A3H (since this product only available after the 2003, as entire year)
    df1a <- df1[df1$date >= "2003-01-01" & df1$date <= "2012-12-31",c("date","modisvar_filled")]
    df1b <- df1a %>% mutate(ymonth = month(date),
                            yday = day(date)) %>% 
      group_by(ymonth, yday) %>% 
      summarise(fpar = mean(modisvar_filled, na.rm = TRUE))
    df1b <- as.data.frame(df1b)[,3] # averaged fapar from 2003 - 2012 (365 length of data)
    df2 <- rep(df1b,(yr_end- yr_start+1)) # repeated it to multiple years, the number of years is consitent to what we collect climate forcing
  } else { # if measurement year bewteen 2003 and 2015 --> use such years directly where consistent with climate forcing
    df1a <- subset(df1,Year>=yr_start & Year<=yr_end)
    df2 <- df1a$modisvar_filled }
  assign(substr(sub('.*daily_', '', fapar_df[i]),1,nchar(sub('.*daily_', '', fapar_df[i]))-4), df2) 
}

#2. Input climate forcing and cbind with fapar
forcing_df <- list.files("/Users/yunpeng/data/NPP_Grassland_final/reprocessing_climate/",full.names = T)
length(forcing_df) # siteinfo_final[375,] doesn't work, since its lat/lon has no measurement data - in edges?
siteinfo_final[375,] 
#forcing - input

for (i in 1:length(forcing_df)){
  df1 <- read.csv(forcing_df[i])
  df1$date <- as.Date(df1$date)
  fapar <- (eval(parse(text=df1$sitename[1])))
  df2 <- cbind(df1,fapar)
  df3 <- df2[,c("date","temp","prec","rain","snow","vpd","ppfd","patm","ccov_int","ccov","fapar","co2")]
  names(df3)[names(df3) == 'rain'] <- 'rainf'
  names(df3)[names(df3) == 'snow'] <- 'snowf'
  assign(paste("final",df1$sitename[1],sep="_"), as_tibble(df3))
}

##extract number of file name: a <- (as.numeric(regmatches(forcing_df[i], regexpr( "\\d+", forcing_df[i]))))
## check if a list of number has lacked a value setdiff(1:485, aa)

#Now, start predicting gpp


#now, we can work on predicting gpp
df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286)


#set all siteinfo_final$whc = 170
siteinfo_final$whc = 170

#now, run for loop for collecting each site - by pass points above that do not have fapar
dim(na_data) #this point needs to be passed
siteinfo_final$pred_gpp_c3 <- NA
siteinfo_final$pred_gpp_c4 <- NA
siteinfo_final$max_vcmax25_c3 <- NA
siteinfo_final$max_vcmax25_c4 <- NA

#now, gpp-c3, gpp-c4, vcmax25max-c3, vcmax25-c4
for (i in 1:nrow(siteinfo_final)) {
  tryCatch({
    #c3 gpp
    forcing <- (eval(parse(text=(paste("final",siteinfo_final$sitename[i],sep="_")))))
    modlist <- run_pmodel_f_bysite( 
      siteinfo_final$sitename[i], 
      params_siml <- list(
        spinup             = TRUE,
        spinupyears        = 10,
        recycle            = 1,
        soilmstress        = TRUE,
        tempstress         = TRUE,
        calc_aet_fapar_vpd = FALSE,
        in_ppfd            = TRUE,
        in_netrad          = FALSE,
        outdt              = 1,
        ltre               = FALSE,
        ltne               = FALSE,
        ltrd               = FALSE,
        ltnd               = FALSE,
        lgr3               = TRUE,
        lgn3               = FALSE,
        lgr4               = FALSE,
        firstyeartrend = siteinfo_final$year_start[i],
        nyeartrend = siteinfo_final$year_end[i]-siteinfo_final$year_start[i]+1), 
      siteinfo = siteinfo_final[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    
    pred_gpp_list <- modlist %>% mutate(ymonth = month(date),yday = day(date)) %>% group_by(ymonth, yday) %>% summarise(gpp = mean(gpp, na.rm = TRUE))
    sum(pred_gpp_list$gpp)
    siteinfo_final[i,c("pred_gpp_c3")] <- sum(pred_gpp_list$gpp)
    
    #c4 gpp
    modlist <- run_pmodel_f_bysite( 
      siteinfo_final$sitename[i], 
      params_siml <- list(
        spinup             = TRUE,
        spinupyears        = 10,
        recycle            = 1,
        soilmstress        = TRUE,
        tempstress         = TRUE,
        calc_aet_fapar_vpd = FALSE,
        in_ppfd            = TRUE,
        in_netrad          = FALSE,
        outdt              = 1,
        ltre               = FALSE,
        ltne               = FALSE,
        ltrd               = FALSE,
        ltnd               = FALSE,
        lgr3               = FALSE,
        lgn3               = FALSE,
        lgr4               = TRUE,
        firstyeartrend = siteinfo_final$year_start[i],
        nyeartrend = siteinfo_final$year_end[i]-siteinfo_final$year_start[i]+1), 
      siteinfo = siteinfo_final[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    pred_gpp_list <- modlist %>% mutate(ymonth = month(date),yday = day(date)) %>% group_by(ymonth, yday) %>% summarise(gpp = mean(gpp, na.rm = TRUE))
    
    siteinfo_final[i,c("pred_gpp_c4")] <- sum(pred_gpp_list$gpp)
    
    #max_vcmax25 - c3
    modlist <- run_pmodel_f_bysite( 
      siteinfo_final$sitename[i], 
      params_siml <- list(
        spinup             = TRUE,
        spinupyears        = 10,
        recycle            = 1,
        soilmstress        = TRUE,
        tempstress         = TRUE,
        calc_aet_fapar_vpd = FALSE,
        in_ppfd            = TRUE,
        in_netrad          = FALSE,
        outdt              = 1,
        ltre               = FALSE,
        ltne               = FALSE,
        ltrd               = FALSE,
        ltnd               = FALSE,
        lgr3               = TRUE,
        lgn3               = FALSE,
        lgr4               = FALSE,
        firstyeartrend = siteinfo_final$year_start[i],
        nyeartrend = siteinfo_final$year_end[i]-siteinfo_final$year_start[i]+1), 
      siteinfo = siteinfo_final[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    max_vcmax25 <- max(modlist$vcmax25)*1000000
    siteinfo_final[i,c("max_vcmax25_c3")] <- max_vcmax25
    
    #max_vcmax25 - c4
    modlist <- run_pmodel_f_bysite( 
      siteinfo_final$sitename[i], 
      params_siml <- list(
        spinup             = TRUE,
        spinupyears        = 10,
        recycle            = 1,
        soilmstress        = TRUE,
        tempstress         = TRUE,
        calc_aet_fapar_vpd = FALSE,
        in_ppfd            = TRUE,
        in_netrad          = FALSE,
        outdt              = 1,
        ltre               = FALSE,
        ltne               = FALSE,
        ltrd               = FALSE,
        ltnd               = FALSE,
        lgr3               = FALSE,
        lgn3               = FALSE,
        lgr4               = TRUE,
        firstyeartrend = siteinfo_final$year_start[i],
        nyeartrend = siteinfo_final$year_end[i]-siteinfo_final$year_start[i]+1), 
      siteinfo = siteinfo_final[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    max_vcmax25 <- max(modlist$vcmax25)*1000000
    siteinfo_final[i,c("max_vcmax25_c4")] <- max_vcmax25
  }, error=function(e){})} 



#collect gpp and combine it into NPP_grassland
siteinfo_final_gpp <- siteinfo_final[,c("sitename","pred_gpp_c3","pred_gpp_c4","max_vcmax25_c3","max_vcmax25_c4")]
names(siteinfo_final_gpp) <- c("sitename2","pred_gpp_c3","pred_gpp_c4","max_vcmax25_c3","max_vcmax25_c4")
NPP_grassland_all2 <-Reduce(function(x,y) merge(x = x, y = y, by = c("sitename2"),all.x=TRUE), 
                            list(NPP_grassland_all,siteinfo_final_gpp))

NPP_grassland_all2 <- NPP_grassland_all2[order(NPP_grassland_all2$no), ]

NPP_grassland$sitename2 <- NPP_grassland_all2$sitename2
NPP_grassland$pred_gpp_c3 <- NPP_grassland_all2$pred_gpp_c3
NPP_grassland$pred_gpp_c4 <- NPP_grassland_all2$pred_gpp_c4
NPP_grassland$max_vcmax25_c3 <- NPP_grassland_all2$max_vcmax25_c3
NPP_grassland$max_vcmax25_c4 <- NPP_grassland_all2$max_vcmax25_c4

summary(siteinfo_final_gpp) # 10 sites missing as expected due to no fapar data available.
subset(siteinfo_final_gpp,is.na(pred_gpp_c3)==TRUE)$sitename2

#input nc file
elev_nc <- read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc")
elev <- as.data.frame(nc_to_df(elev_nc, varnam = "elevation"))

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

#cbind all predictors, and its lon, lat, z
all_predictors <- cbind(elev,Tg$myvar,PPFD$myvar,vpd$myvar,
                        alpha$myvar,fAPAR$myvar,age$myvar)

names(all_predictors) <- c("lon","lat","z","Tg","PPFD","vpd",
                           "alpha","fAPAR","age")

Tg_df <- all_predictors[,c("lon","lat","z","Tg")]
PPFD_df <- all_predictors[,c("lon","lat","z","PPFD")]
vpd_df <- all_predictors[,c("lon","lat","z","vpd")]
alpha_df <- all_predictors[,c("lon","lat","z","alpha")]
fAPAR_df <- all_predictors[,c("lon","lat","z","fAPAR")]
age_df <- all_predictors[,c("lon","lat","z","age")]

#now, apply gwr to extract site predictors' value
grassland_site <- NPP_grassland[,c("lon","lat","z")]
grassland_site$Tg <- NA
grassland_site$PPFD <- NA
grassland_site$vpd <- NA
grassland_site$alpha <- NA
#grassland_site$age <- NA
grassland_site$fapar <- NA

a <- 1.5 # which degree (distance) of grid when interpolating gwr from global grids

#Extract Tg, PPFD, vpd, alpha,fAPAR

for (i in 1:nrow(grassland_site)) {
  tryCatch({
    #Tg
    Tg_global <- na.omit(Tg_df)
    NRE_part <- subset(Tg_global,lon>(grassland_site[i,1]-a)&lon<(grassland_site[i,1]+a)&
                         lat>(grassland_site[i,2]-a)&lat<(grassland_site[i,2]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- grassland_site[i,1:3]
    coordinates(NRE_coord) <- c("lon","lat")
    grassland_site[i,c("Tg")] <- (gwr(Tg ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #ppfd
    PPFD_global <- na.omit(PPFD_df)
    NRE_part <- subset(PPFD_global,lon>(grassland_site[i,1]-a)&lon<(grassland_site[i,1]+a)&
                         lat>(grassland_site[i,2]-a)&lat<(grassland_site[i,2]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- grassland_site[i,1:3]
    coordinates(NRE_coord) <- c("lon","lat")
    grassland_site[i,c("PPFD")] <- (gwr(PPFD ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #vpd
    vpd_global <- na.omit(vpd_df)
    NRE_part <- subset(vpd_global,lon>(grassland_site[i,1]-a)&lon<(grassland_site[i,1]+a)&
                         lat>(grassland_site[i,2]-a)&lat<(grassland_site[i,2]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- grassland_site[i,1:3]
    coordinates(NRE_coord) <- c("lon","lat")
    grassland_site[i,c("vpd")] <- (gwr(vpd ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #alpha
    alpha_global <- na.omit(alpha_df)
    NRE_part <- subset(alpha_global,lon>(grassland_site[i,1]-a)&lon<(grassland_site[i,1]+a)&
                         lat>(grassland_site[i,2]-a)&lat<(grassland_site[i,2]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- grassland_site[i,1:3]
    coordinates(NRE_coord) <- c("lon","lat")
    grassland_site[i,c("alpha")]  <- (gwr(alpha ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #fAPAR
    fAPAR_global <- na.omit(fAPAR_df)
    NRE_part <- subset(fAPAR_global,lon>(grassland_site[i,1]-a)&lon<(grassland_site[i,1]+a)&
                         lat>(grassland_site[i,2]-a)&lat<(grassland_site[i,2]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- grassland_site[i,1:3]
    coordinates(NRE_coord) <- c("lon","lat")
    grassland_site[i,c("fapar")]<- (gwr(fAPAR ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
  }, error=function(e){})} 


summary(grassland_site)

library(raster)
library(rgdal)
library(dplyr)
library(rbeni)
library(ncdf4)
soil <- raster('~/data/ISRIC/data_orig/data/raster/w001000.adf')
NRE_lonlat <- grassland_site[,c("lon","lat","z")]

sp_sites <- SpatialPoints(NRE_lonlat[,c("lon","lat","z")]) # only select lon and lat

#change its variable name to SUID, this is a unique code that could be used to merged with soil data, which will be further merged with csv below.
NRE_lonlat2 <- raster::extract(soil, sp_sites, sp = TRUE) %>% as_tibble() %>% 
  right_join(NRE_lonlat, by = c("lon", "lat","z")) %>% 
  dplyr::rename( SUID = w001000)

#input soil information data csv
ISRIC.data<-read.csv(file="~/data/ISRIC/data_orig/data/HW30s_FULL.csv",header=TRUE,sep=";",dec = ".") # Now, input ISRIC database
data.soil.extract <- merge(NRE_lonlat2,ISRIC.data,by='SUID',all.x=TRUE) # merge site with soil variables by using SUID
data.soil.extract2 <- subset(data.soil.extract,CNrt>0) # select available CNrt
data.soil.extract3 <- data.soil.extract2[,c("lon","lat","z","CNrt")] # only select CNrt variable

# note that in each site there might be more than 1 samples measured, so we should aggregate them which make sures that one grid holds one data only.
ss1 <- aggregate(data.soil.extract3,by=list(data.soil.extract3$lon,data.soil.extract3$lat,data.soil.extract3$z), FUN=mean, na.rm=TRUE) 
ss2 <- ss1[,c("lon","lat","z","CNrt")] # now,select lon, lat, z and CNrt only

# finally, merging site-based soil c/n data into our current dataframe
grassland_site$no <- c(1:nrow(grassland_site))
grassland_site2 <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat","z"),all.x=TRUE), 
                         list(grassland_site,ss2))

grassland_site3 <- grassland_site2[order(grassland_site2$no), ]

head(grassland_site3)

NPP_grassland_final <- cbind(NPP_grassland,grassland_site3[,c(4,5,6,7,8,10)])
summary(NPP_grassland_final)

#now, process final data, and used in final simulations
NPP_grassland_final$fapar[NPP_grassland_final$fapar>1] <- NA
NPP_grassland_final$alpha[NPP_grassland_final$alpha>1] <- NA
#NPP_grassland_final$age[NPP_grassland_final$age<0] <- NA
summary(NPP_grassland_final)



##########NOW, final calculating gpp based on weighted sum method --> the first is only based on measured c3c4, the second is the combination of measured + map c3c4.
NPP_grassland_final$weightedgpp_measured_c3 <- (NPP_grassland_final$pred_gpp_c3 * NPP_grassland_final$c3_percentage)+(NPP_grassland_final$pred_gpp_c4 * (1-NPP_grassland_final$c3_percentage))
NPP_grassland_final$weightedgpp_all <- (NPP_grassland_final$pred_gpp_c3 * NPP_grassland_final$c3_percentage_final)+(NPP_grassland_final$pred_gpp_c4 * (1-NPP_grassland_final$c3_percentage_final))

#firstly, gpp
#replace 3 sites predicted GPP to stocker et al. 2020 GMD GPP (where using fluxnet measured climate data)
#the spatial based gpp prediction was available in /Users/yunpeng/data/gpp_gmd/gpp_pmodel_fluxnet2015_stocker19gmd_spatial.csv
#...and also in zenedo: https://zenodo.org/record/3559850#.YHgs4pMzaqA
#Stocker et al. 2019 GMD: Difference between ORG (Wang et al. 2017), BRC (newly have temperature-dependent phio) and FULL (BRC + soil moisture correction).
#our rsofun is FULL (T-dependent phio + soil moisture correction) so we replace it also to FULL version data
#our current simulation
NPP_grassland_final$weightedgpp_all[NPP_grassland_final$site=="CN-du2-D01"]
NPP_grassland_final$weightedgpp_all[NPP_grassland_final$site=="RU-ha1-F01"]
NPP_grassland_final$weightedgpp_all[NPP_grassland_final$site=="DE-gri-D01"]
#change to (see csv location - site-name at FULL setup - not changing too much!)
NPP_grassland_final$weightedgpp_all[NPP_grassland_final$site=="CN-du2-D01"] <- 541.8902364
NPP_grassland_final$weightedgpp_all[NPP_grassland_final$site=="RU-ha1-F01"] <- 536.0214562
NPP_grassland_final$weightedgpp_all[NPP_grassland_final$site=="DE-gri-D01"] <- 1491.184888

#convert management code
NPP_grassland_final$Management.code[is.na(NPP_grassland_final$Management.code)==TRUE] <- "N"

#now, summarise by pft, management and file - and filter.
NPP_grassland_final %>% group_by(pft)  %>% summarise(number = n())
NPP_grassland_final2 <- subset(NPP_grassland_final,pft=="Grassland")

NPP_grassland_final2 %>% group_by(Management.code)  %>% summarise(number = n())
NPP_grassland_final3 <- subset(NPP_grassland_final2,Management.code=="N" | Management.code== "SN")

# make Tiandi's bnpp and npp = NA
NPP_grassland_final3 %>% group_by(file)  %>% summarise(number = n())

NPP_grassland_final3$BNPP_1[NPP_grassland_final3$file=="Tiandi Grassland"] <- NA
NPP_grassland_final3$TNPP_1[NPP_grassland_final3$file=="Tiandi Grassland"] <- NA

NPP_grassland_final4 <- NPP_grassland_final3[,!(names(NPP_grassland_final3) %in% c("NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","Management"))]

subset(NPP_grassland_final4,GPP>0)
#Now, time to examine our data
My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size = 20))


NPP_grassland_final5_gpp_npp_anpp <- aggregate(NPP_grassland_final4,by=list(NPP_grassland_final4$lon,NPP_grassland_final4$lat,NPP_grassland_final4$z), FUN=mean, na.rm=TRUE)

#npp/gpp model. Now, sites are too less when construct model for npp/gpp, let's aggregate them firstly.
tnpp_grass <- (lm((TNPP_1)~-1+(weightedgpp_all),data=NPP_grassland_final5_gpp_npp_anpp)) #0.435 for using weighted gpp (df = 79)
summary(tnpp_grass)
#summary(lm((TNPP_1)~-1+(GPP),data=NPP_grassland_final5_gpp_npp_anpp)) #0.389 for using measured gpp (df=18)
save(tnpp_grass, file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/tnpp_grass.RData")

#anpp/gpp model.
anpp_grass <- (lm((ANPP_2)~-1+(weightedgpp_all),data=NPP_grassland_final5_gpp_npp_anpp)) #anpp/tnpp = 0.44730 (df = 78)
summary(anpp_grass)
save(anpp_grass, file = "/Users/yunpeng/data/NPP_grassland_final/statistical_model/anpp_grass.RData")

#leaf c/n model. median = 18
summary(NPP_grassland_final5_gpp_npp_anpp$CN_leaf_final)

#root c/n model median = 41
summary(NPP_grassland_final5_gpp_npp_anpp$CN_root_final)

#constant only
NPP_grassland_final4$pred_npp <- NPP_grassland_final4$weightedgpp_all * 0.435
NPP_grassland_final4$pred_anpp <- NPP_grassland_final4$weightedgpp_all * 0.228
NPP_grassland_final4$pred_lnf <- NPP_grassland_final4$pred_anpp/18
NPP_grassland_final4$pred_bnpp <- NPP_grassland_final4$pred_npp - NPP_grassland_final4$pred_anpp
NPP_grassland_final4$pred_bnf <- NPP_grassland_final4$pred_bnpp/41

#anpp_gpp <- 0.302399 + 0.0126378*NPP_grassland_final4$Tg + -0.580933*NPP_grassland_final4$fapar
#anpp_gpp[anpp_gpp<0]
#anpp_gpp[anpp_gpp<0] <- NA
#NPP_grassland_final4$pred_anpp_reg <- anpp_gpp * NPP_grassland_final4$weightedgpp_all

#gpp r2 = 0.29
NPP_grassland_final4
csvfile <- paste("/Users/yunpeng/data/NPP_Grassland_final/NPP_grass_validation.csv")
write.csv(NPP_grassland_final4, csvfile, row.names = TRUE)

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 20))

ggplot(NPP_grassland_final4, aes(x=weightedgpp_all, y=GPP)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic() + My_Theme #+ggtitle("Observed GPP vs. Predicted GPP")
summary(lm(GPP~weightedgpp_all,data=NPP_grassland_final4))

#####npp r2 = 0.1128
ggplot(NPP_grassland_final4, aes(x=pred_npp, y=TNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic() + My_Theme 
summary(lm(TNPP_1~pred_npp,data=NPP_grassland_final4))

####anpp r2 = 0.1732
ggplot(NPP_grassland_final4, aes(x=pred_anpp, y=ANPP_2)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic() + My_Theme 
summary(lm(ANPP_2~pred_anpp,data=NPP_grassland_final4)) 

#leaf N flux r2 = 0.096
ggplot(NPP_grassland_final4, aes(x=pred_lnf, y=lnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+ xlim(c(0,10))+
  xlab("Prediction")+ylab("Observation")+theme_classic() + My_Theme
summary(lm(lnf_obs_final~pred_lnf,data=NPP_grassland_final4)) 

#add legume and biome
legume_N <- read.csv("/Users/yunpeng/data/npp_stoichiometry_grasslands_tiandi/China_grassland_CN_stoichiometry_with_matched_NPP_species_legume_20201214.csv")
legume_N_only <- subset(legume_N,Legume_CN=="N"|Legume_CN=="N_N"|
                          Legume_CN=="N_N_N"|Legume_CN=="N_N_N_N"|Legume_CN=="N_N_N_N_N")
legume_N_only <- legume_N_only[,c("Longitude_CN","Latitude_CN","Altitude_CN")]
legume_N_only_site <- aggregate(legume_N_only,by=list(legume_N_only$Longitude_CN,
                                                      legume_N_only$Latitude_CN,
                                                      legume_N_only$Altitude_CN), mean,na.rm=TRUE)
legume_N_only_site <- legume_N_only_site[,c("Longitude_CN","Latitude_CN","Altitude_CN")]
names(legume_N_only_site) <- c("lon","lat","z")
legume_N_only_site$legume <- "N"
NPP_grassland_final4$no <- c(1:nrow(NPP_grassland_final4))
NPP_grassland_final4_legume <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat","z"),all.x=TRUE), 
                                     list(NPP_grassland_final4,legume_N_only_site))
NPP_grassland_final4_legume <- NPP_grassland_final4_legume[order(NPP_grassland_final4_legume$no), ]

subset(NPP_grassland_final4_legume,lnf_obs_final>0) %>% group_by(legume) %>% summarise(number = n())
summary(subset(NPP_grassland_final4_legume,legume=="N")$CN_leaf_final)
summary(subset(NPP_grassland_final4_legume,is.na(legume)==TRUE)$CN_leaf_final)

ggplot(NPP_grassland_final4_legume, aes(x=pred_lnf, y=lnf_obs_final)) +
  geom_point(aes(color=factor(legume)))+geom_abline(intercept=0,slope=1)+geom_smooth(aes(color=factor(legume)),method = "lm", se = TRUE)+ xlim(c(0,10))+
  xlab("Predicted leaf N flux")+ylab("Measured leaf N flux")+theme_classic() + My_Theme
summary(lm(lnf_obs_final~pred_lnf,data=NPP_grassland_final4)) 
summary(lm(lnf_obs_final~pred_lnf,data=subset(NPP_grassland_final4_legume,legume=="N"))) 

china <- (subset(NPP_grassland_final4_legume,file=="Tiandi Grassland"))

china <- china[order(china$lon,china$lat,china$z,china$ANPP_2,china$CN_leaf_final), ]
dim(china)

biome_more <- read.csv("/Users/yunpeng/data/npp_stoichiometry_grasslands_tiandi/npp_stoichiometry_china_grassland_CN_stoichiometry_with_matched_NPP_data_from_Prof_Fang_group_20201026.csv")
dim(biome_more)
biome_more <- biome_more[order(biome_more$Longitude_stoichiometry,biome_more$Latitude_stoichiometry,
                               biome_more$Altitude_stoichiometry,biome_more$ANPP,biome_more$CNratio_leaf), ]
summary(biome_more$Longitude_stoichiometry - china$lon)
summary(biome_more$Latitude_stoichiometry - china$lat)
summary(biome_more$Altitude_stoichiometry - china$z)

china$Vegetation_type_stoichiometry <- biome_more$Vegetation_type_stoichiometry

china$pred_lnf_new <- china$weightedgpp_all*0.228/18
ggplot(china, aes(x=pred_lnf, y=lnf_obs_final)) +geom_abline(intercept=0,slope=1)+geom_smooth(aes(color=factor(Vegetation_type_stoichiometry)),method = "lm", se = TRUE)+ xlim(c(0,10))+
  xlab("Predicted leaf N flux")+ylab("Measured leaf N flux")+theme_classic() + My_Theme

#do some filtering here
china %>% group_by(Vegetation_type_stoichiometry)  %>% summarise(number = n())
china$pft_china[china$Vegetation_type_stoichiometry=="Alpine meadow"|
                  china$Vegetation_type_stoichiometry=="Alpine scrub"|
                  china$Vegetation_type_stoichiometry=="Alpine steppe"] <- "Alpine (meadow, scrub, steppe) biome"

china$pft_china[china$Vegetation_type_stoichiometry=="Desert"|
                  china$Vegetation_type_stoichiometry=="Desert steppe"] <- "Desert biome"

china$pft_china[china$Vegetation_type_stoichiometry=="Meadow steppe"|
                  china$Vegetation_type_stoichiometry=="Typical steppe"] <- "Typical and meadow steppe"  

china$pft_china[china$Vegetation_type_stoichiometry=="Temperate meadow"|
                  china$Vegetation_type_stoichiometry=="Temperate scrub"] <- "Temperate (meadow, scrub) biome"  

china %>% group_by(c3_percentage_final)  %>% summarise(number = n())

china$c3_china[china$c3_percentage_final<1] <- "c4"
china$c3_china[china$c3_percentage_final==1] <- "c3"

china$c3_china[china$c3_percentage_final<0.5] <- "c4"
china$c3_china[china$c3_percentage_final>=0.5] <- "c3"


ggplot(china, aes(x=pred_lnf, y=lnf_obs_final)) +geom_abline(intercept=0,slope=1)+geom_smooth(aes(color=factor(c3_percentage_final)),method = "lm", se = TRUE)+ xlim(c(0,10))+
  xlab("Predicted leaf N flux")+ylab("Measured leaf N flux")+theme_classic() + My_Theme

#first - c3 and c4
#using plyr to produce means for each type
library(plyr)
ggplot(china, aes(x = lnf_obs_final, fill = c3_china)) +labs(x = "Leaf N uptake (gN/m2/yr)")+
  geom_density(alpha = .3) + #alpha used for filling the density
  geom_vline(data =  ddply(china, "c3_china", summarise, rating.mean = mean(lnf_obs_final,na.rm=TRUE)),
             aes(xintercept = rating.mean, colour = c3_china),
             linetype = "longdash", size=1)+ theme_classic()+ My_Theme+ theme(legend.title = element_text(size = 20),legend.text = element_text(size = 15))

ggplot(china, aes(x = lnf_obs_final, fill = pft_china)) +labs(x = "Leaf N uptake (gN/m2/yr)")+
  geom_density(alpha = .3) + #alpha used for filling the density
  geom_vline(data = ddply(china, "pft_china", summarise, rating.mean = mean(lnf_obs_final,na.rm=TRUE)), 
             aes(xintercept = rating.mean, colour = pft_china),
             linetype = "longdash", size=1)+ theme_classic()+ My_Theme+ theme(legend.title = element_text(size = 20),legend.text = element_text(size = 15))

china$legume[china$legume=="N"] <- "Non-legume species"
china$legume[is.na(china$legume)==TRUE] <- "Legume species"

ggplot(china, aes(x = lnf_obs_final, fill = legume)) +labs(x = "Leaf N uptake (gN/m2/yr)")+
  geom_density(alpha = .3) + #alpha used for filling the density
  geom_vline(data = ddply(china, "legume", summarise, rating.mean = mean(lnf_obs_final,na.rm=TRUE)), 
             aes(xintercept = rating.mean, colour = legume),
             linetype = "longdash", size=1)+ theme_classic()+ My_Theme+ theme(legend.title = element_text(size = 20),legend.text = element_text(size = 15))

#bnpp

ggplot(NPP_grassland_final4, aes(x=pred_bnpp, y=BNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+ 
  xlab("Predicted BNPP")+ylab("Measured BNPP")+theme_classic() + My_Theme
summary(lm(BNPP_1~pred_bnpp,data=NPP_grassland_final4)) 

#bnf
ggplot(NPP_grassland_final4, aes(x=pred_bnf, y=bnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+ 
  xlab("Predicted BNF")+ylab("Measured BNF")+theme_classic() + My_Theme
summary(lm(bnf_obs_final~pred_bnf,data=NPP_grassland_final4)) 

save.image(file = "/Users/yunpeng/data/NPP_Grassland_final/grass_simulation.Rdata")
