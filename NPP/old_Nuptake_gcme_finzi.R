rm(list=ls())
#library(ingestr)
library(dplyr)
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
library(spgwr)
library(maps)
library(rworldmap)
library(cowplot)
devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/rsofun/")

#load(file = "~/yunkepeng/nimpl_sofun_inputs/forest/New_Nuptake_site_simulation.Rdata")

#### Input N uptake
#(1) newly added Nmin rate data from Finzi
Finzi <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nmin_Finzi/Nmin_Finzi.csv")
names(Finzi)[names(Finzi) == "Lat"] <- "lat"
names(Finzi)[names(Finzi) == "Long"] <- "lon"
devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/ingestr/")

#See Gill and Finzi Fig3 for pft info 
#Grassland
#Finzi_Grassland <- subset(Finzi, Biome=="temp grass")
#Finzi_Grassland_sitemean <- aggregate(Finzi_Grassland,by=list(Finzi_Grassland$lon,Finzi_Grassland$lat), FUN=mean, na.rm=TRUE) #site-mean
#dim(Finzi_Grassland_sitemean)
#for (i in 1:nrow(Finzi_Grassland_sitemean)){
#  Finzi_Grassland_sitemean$sitename[i] <- paste("Finzi_Grass",i,sep = "") # this is also sitename for fpar
#  Finzi_Grassland_sitemean$sitename_climate[i] <- paste("Finzi_Grass_climate",i,sep = "")
#}
#df_etopo <- ingest(Finzi_Grassland_sitemean,source = "etopo1",dir = "~/data/etopo/" )
#Finzi_Grassland_sitemean$elv <- as.numeric(as.data.frame(df_etopo$data))
#Finzi_Grassland_sitemean

#Forest - only merging forest this time
Finzi_Forest <- subset(Finzi, Biome!="temp grass")
Finzi_Forest_sitemean <- aggregate(Finzi_Forest,by=list(Finzi_Forest$lon,Finzi_Forest$lat), FUN=mean, na.rm=TRUE) #site-mean
dim(Finzi_Forest_sitemean)
for (i in 1:nrow(Finzi_Forest_sitemean)){
  Finzi_Forest_sitemean$sitename[i] <- paste("Finzi_Forest",i,sep = "") # this is also sitename for fpar
  Finzi_Forest_sitemean$sitename_climate[i] <- paste("Finzi_Forest_climate",i,sep = "")
  
}
df_etopo <- ingest(Finzi_Forest_sitemean,source = "etopo1",dir = "~/data/etopo/" )
Finzi_Forest_sitemean$elv <- as.numeric(as.data.frame(df_etopo$data))
Finzi_Forest_sitemean$elv[Finzi_Forest_sitemean$elv< 0] <- 0
Finzi_Forest_sitemean
#Finzi_final <- dplyr::bind_rows(Finzi_Forest_sitemean, Finzi_Grassland_sitemean)

Finzi_Forest_sitemean2 <- Finzi_Forest_sitemean[,c("lon","lat","elv","sitename","sitename_climate")]
dim(Finzi_Forest_sitemean2)
Finzi_all <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat"),all.x=TRUE), 
                   list(Finzi,Finzi_Forest_sitemean2))

Finzi_all_forest <- subset(Finzi_all, Biome!="temp grass" & is.na(lon)==FALSE)
summary(Finzi_all_forest)
Finzi_all_forest$year_start <- 1984
Finzi_all_forest$year_end <- 2013

#(2) Nuptake from gcme
gcme_nuptake <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nuptake_gcme/gcme_nuptake_coord_interpolated.csv")
for (i in 1:nrow(gcme_nuptake)){
  gcme_nuptake$sitename[i] <- paste("gcme",i,sep = "") # this is also sitename for fpar
  gcme_nuptake$sitename_climate[i] <- paste("gcme_climate",i,sep = "")
}

#elv
df_etopo <- ingest(gcme_nuptake,source = "etopo1",dir = "~/data/etopo/" )
gcme_nuptake$elv <- as.numeric(as.data.frame(df_etopo$data))
gcme_nuptake$elv[gcme_nuptake$elv<0] <- 0

siteinfo_gcme <- data.frame(
  sitename = gcme_nuptake$sitename,
  sitename_climate = gcme_nuptake$sitename_climate,
  lon = gcme_nuptake$lon,
  lat = gcme_nuptake$lat,
  elv = gcme_nuptake$elv,
  year_start = gcme_nuptake$start_yr,
  year_end = gcme_nuptake$start_yr + gcme_nuptake$Year_long -1
)

siteinfo_gcme$exp_nam <- gcme_nuptake$exp_nam

gcme_data <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nuptake_gcme/gcme_nuptake_data.csv")
gcme_data <- gcme_data[,c("exp_nam","factors","treatment","ambient","elevated","Unit","id","Start_Year","Year","exp_nam")]

#all converting to gN/m2/yr
gcme_data$ambient[gcme_data$Unit=="Kg_N_ha-1"] <- gcme_data$ambient[gcme_data$Unit=="Kg_N_ha-1"]/10
gcme_data$ambient[gcme_data$Unit=="kg_N/ha"] <- gcme_data$ambient[gcme_data$Unit=="kg_N/ha"]/10
gcme_data$ambient[gcme_data$Unit=="mg_N/kg*day_"] <- NA #quite weired about the unit, for one site. Disregard them first
hist(gcme_data$ambient)

gcme_data$elevated[gcme_data$Unit=="Kg_N_ha-1"] <- gcme_data$elevated[gcme_data$Unit=="Kg_N_ha-1"]/10
gcme_data$elevated[gcme_data$Unit=="kg_N/ha"] <- gcme_data$elevated[gcme_data$Unit=="kg_N/ha"]/10
gcme_data$elevated[gcme_data$Unit=="mg_N/kg*day_"] <- NA #quite weired about the unit, for one site. Disregard them first
hist(gcme_data$elevated)

#why some values are too low? - but we don't remove them now
gcme_data <- subset(gcme_data,ambient>0)
subset(gcme_data,ambient<1) %>% group_by(exp_nam) %>% summarise(number = n())
#after look, "RiceFACE_Japan_A_1998_39,40_141" looks fine, as most of them were still in good range. But remove the other 5 sites (they may be wrong due to measurement mistakes or unit errors, we don't know)
#gcme_data_final <- subset(gcme_data,exp_nam!="Michigan_UNDERC_bog" & exp_nam!="Michigan_UNDERC_intermFen" & exp_nam!="Michigan_UNDERC_richFen"& exp_nam!="RiceFACE_China_32N_120E_Or_Tr_7"& exp_nam!="TL_7")
gcme_data_final <- gcme_data
hist(gcme_data_final$ambient)

#finally, mergeing them to obtain siteinfo

gcme_data_final_forest <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam"),all.x=TRUE),list(gcme_data_final,siteinfo_gcme))

#rbind them to get final Nmin. data

Finzi_all_forest <- Finzi_all_forest[,c("lon","lat","elv","sitename","sitename_climate","year_start","year_end","Nmin")]
Finzi_all_forest <- rename(Finzi_all_forest, obs_nuptake = Nmin)
Finzi_all_forest$method <- "Net minerlization (Finzi paper)"
gcme_data_final_forest <- gcme_data_final_forest[,c("lon","lat","elv","sitename","sitename_climate","year_start","year_end","ambient",
                                                    "elevated","factors","treatment","Unit","id","Start_Year","Year","exp_nam")]
gcme_data_final_forest <- rename(gcme_data_final_forest, obs_nuptake = ambient)
gcme_data_final_forest$method <- "Total Nuptake reported in GCME, at ambient condition"

Nmin_final <- dplyr::bind_rows(Finzi_all_forest,gcme_data_final_forest)
Nmin_final

aa <- subset(Nmin_final,is.na(exp_nam)==FALSE) %>% group_by(exp_nam) %>% summarise (number=n())

#input Filzi's climates
forcing_df <- list.files("~/data/NPP_Yunke/Nmin_Finzi/reprocessing_Nmin/climates/",full.names = T) # 2 points were missing, as expected
length(forcing_df)

fapar_df <- list.files("~/data/NPP_Yunke/Nmin_Finzi/reprocessing_Nmin/fapar/",full.names = T)
length(fapar_df)-1

fapar_org_df <- list.files("~/data/NPP_Yunke/Nmin_Finzi/reprocessing_Nmin/fapar/raw/",full.names = T)
length(fapar_org_df)

#1. fapar - input

#1. fpar - check missing data - and also, input all years fapar (2001-2015), which will be selected in measurement year only later on 
for (i in 1:(length(fapar_df)-1)){
  df1 <- read.csv(fapar_df[i])
  df1$date <- as.Date(df1$date)
  df1 <- df1[!(format(df1$date,"%m") == "02" & format(df1$date, "%d") == "29"), , drop = FALSE]
  df2 <- df1[,c("date","modisvar_filled")]
  assign(substr(sub('.*daily_', '', fapar_df[i]),1,nchar(sub('.*daily_', '', fapar_df[i]))-4), df2) 
}

#2. forcing - combing fapar and climates into a df.
for (i in 1:(length(forcing_df)-18)){ # remove the later 18 sites where is grassland
  df1 <- read.csv(forcing_df[i])
  df1$date <- as.Date(df1$date)
  
  sitename_climate <- subset(Nmin_final,Nmin_final$sitename_climate == df1$sitename[1])$sitename_climate
  sitename_fapar <- subset(Nmin_final,Nmin_final$sitename_climate == df1$sitename[1])$sitename
  
  fapar <- (eval(parse(text=sitename_fapar)))
  fapar$Year <- year(fapar$date)
  
  yr_start <- 1984
  yr_end <- 2013
  
  if (yr_start<=2002) { # if measurement year before 2002 for a certain site --> calculating 2003-2012 average of MCD15A3H (since this product only available after the 2003, as entire year)
    df1a <- fapar[fapar$date >= "2003-01-01" & fapar$date <= "2012-12-31",c("date","modisvar_filled")]
    df1b <- df1a %>% mutate(ymonth = month(date),
                            yday = day(date)) %>% 
      group_by(ymonth, yday) %>% 
      summarise(fpar = mean(modisvar_filled, na.rm = TRUE))
    df1b <- as.data.frame(df1b)[,3] # averaged fapar from 2003 - 2012 (365 length of data)
    df2 <- rep(df1b,(yr_end- yr_start+1)) # repeated it to multiple years, the number of years is consitent to what we collect climate forcing
  } else { # if measurement year bewteen 2003 and 2015 --> use such years directly where consistent with climate forcing
    df1a <- subset(fapar,Year>=yr_start & Year<=yr_end)
    df2 <- df1a$modisvar_filled }
  
  fpar <- df2
  
  df3 <- cbind(df1,fpar)
  df3 <- df3[,c("date","temp","prec","rain","snow","vpd","ppfd","patm","ccov_int","ccov","fpar","co2")]
  names(df3)[names(df3) == 'rain'] <- 'rainf'
  names(df3)[names(df3) == 'snow'] <- 'snowf'
  names(df3)[names(df3) == 'fpar'] <- 'fapar'
  assign(paste("final",df1$sitename[1],sep="_"), as_tibble(df3))
}

#3. rsofun to predict gpp
df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286)

Nmin_final$whc = 170

Nmin_final$pred_gpp_c3 <- NA
#NPP_Forest$pred_gpp_c4 <- NA
Nmin_final$max_vcmax25_c3 <- NA
#NPP_Forest$max_vcmax25_c4 <- NA


#using rsofun
for (i in 1:nrow(Nmin_final)) {
  tryCatch({
    #c3
    forcing <- (eval(parse(text=(paste("final",Nmin_final$sitename_climate[i],sep="_")))))
    modlist <- run_pmodel_f_bysite( 
      Nmin_final$sitename_climate[i], 
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
        firstyeartrend = Nmin_final$year_start[i],
        nyeartrend = Nmin_final$year_end[i]-Nmin_final$year_start[i]+1), 
      siteinfo = Nmin_final[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    
    pred_gpp_list <- modlist %>% mutate(ymonth = month(date),yday = day(date)) %>% group_by(ymonth, yday) %>% summarise(gpp = mean(gpp, na.rm = TRUE))
    max_vcmax25 <- max(modlist$vcmax25)*1000000
    
    Nmin_final[i,c("pred_gpp_c3")] <- sum(pred_gpp_list$gpp)
    Nmin_final[i,c("max_vcmax25_c3")] <- max_vcmax25
  }, error=function(e){})} 
#for fapar - if not available for n_focal = 0, then changing to 1, then 2...


#input gcme climates
forcing_df <- list.files("~/data/NPP_Yunke/Nuptake_gcme/climates/",full.names = T) # 2 points were missing, as expected due to 
length(forcing_df)

fapar_df <- list.files("~/data/NPP_Yunke/Nuptake_gcme/fapar/",full.names = T)
length(fapar_df)-1

#1. fapar - input

#1. fpar - check missing data - and also, input all years fapar (2001-2015), which will be selected in measurement year only later on 
for (i in 1:(length(fapar_df)-1)){
  df1 <- read.csv(fapar_df[i])
  df1$date <- as.Date(df1$date)
  df1 <- df1[!(format(df1$date,"%m") == "02" & format(df1$date, "%d") == "29"), , drop = FALSE]
  df2 <- df1[,c("date","modisvar_filled")]
  assign(substr(sub('.*daily_', '', fapar_df[i]),1,nchar(sub('.*daily_', '', fapar_df[i]))-4), df2) 
}

#2. forcing - combing fapar and climates into a df.
for (i in 1:(length(forcing_df))) {
  tryCatch({
    df1 <- read.csv(forcing_df[i])
    df1$date <- as.Date(df1$date)
    
    sitename_climate <- subset(Nmin_final,Nmin_final$sitename_climate == df1$sitename[1])$sitename_climate[1]
    sitename_fapar <- subset(Nmin_final,Nmin_final$sitename_climate == df1$sitename[1])$sitename[1]
    
    fapar <- (eval(parse(text=sitename_fapar)))
    fapar$Year <- year(fapar$date)
    
    yr_start <- subset(Nmin_final,Nmin_final$sitename_climate == df1$sitename[1])$year_start[1]
    yr_end <- subset(Nmin_final,Nmin_final$sitename_climate == df1$sitename[1])$year_end[1]
    
    if (yr_start<=2002) { # if measurement year before 2002 for a certain site --> calculating 2003-2012 average of MCD15A3H (since this product only available after the 2003, as entire year)
      df1a <- fapar[fapar$date >= "2003-01-01" & fapar$date <= "2012-12-31",c("date","modisvar_filled")]
      df1b <- df1a %>% mutate(ymonth = month(date),
                              yday = day(date)) %>% 
        group_by(ymonth, yday) %>% 
        summarise(fpar = mean(modisvar_filled, na.rm = TRUE))
      df1b <- as.data.frame(df1b)[,3] # averaged fapar from 2003 - 2012 (365 length of data)
      df2 <- rep(df1b,(yr_end- yr_start+1)) # repeated it to multiple years, the number of years is consitent to what we collect climate forcing
    } else { # if measurement year bewteen 2003 and 2015 --> use such years directly where consistent with climate forcing
      df1a <- subset(fapar,Year>=yr_start & Year<=yr_end)
      df2 <- df1a$modisvar_filled }
    
    fpar <- df2
    
    df3 <- cbind(df1,fpar)
    df3 <- df3[,c("date","temp","prec","rain","snow","vpd","ppfd","patm","ccov_int","ccov","fpar","co2")]
    names(df3)[names(df3) == 'rain'] <- 'rainf'
    names(df3)[names(df3) == 'snow'] <- 'snowf'
    names(df3)[names(df3) == 'fpar'] <- 'fapar'
    assign(paste("final",df1$sitename[1],sep="_"), as_tibble(df3))
  }, error=function(e){})} 


#3. rsofun to predict gpp
df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286)

Nmin_final$whc = 170

#using rsofun for gcme sites only
for (i in 241:429) {
  tryCatch({
    #c3
    forcing <- (eval(parse(text=(paste("final",Nmin_final$sitename_climate[i],sep="_")))))
    modlist <- run_pmodel_f_bysite( 
      Nmin_final$sitename_climate[i], 
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
        firstyeartrend = Nmin_final$year_start[i],
        nyeartrend = Nmin_final$year_end[i]-Nmin_final$year_start[i]+1), 
      siteinfo = Nmin_final[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    
    pred_gpp_list <- modlist %>% mutate(ymonth = month(date),yday = day(date)) %>% group_by(ymonth, yday) %>% summarise(gpp = mean(gpp, na.rm = TRUE))
    max_vcmax25 <- max(modlist$vcmax25)*1000000
    
    Nmin_final[i,c("pred_gpp_c3")] <- sum(pred_gpp_list$gpp)
    Nmin_final[i,c("max_vcmax25_c3")] <- max_vcmax25
  }, error=function(e){})} 

Nmin_final

NPP_Forest <- Nmin_final
#already checked here -all NA in pred_gpp_c3 is due to NA of fAPAR in orig/, either in nfocal = 0, 1 and 2.

#stop at here and save - where "NPP_Forest_all_flux" can be used furtherly, without the need to re-load

#now, extracting site values from Tg, alpha, c/n.....
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

CNrt <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/CNrt.nc"),
  varnam = "CNrt"))

LMA <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/LMA.nc"),
  varnam = "LMA"))

#2. Create a function to specify path, loop many years nc file and output a dataframe (lon, lat, var).
inputnc <- function(name,start_year,end_year){
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

#cbind all predictors, and its lon, lat, z
all_predictors <- cbind(elev,Tg$myvar,PPFD$myvar,vpd$myvar,
                        alpha$myvar,fAPAR$myvar,age$myvar,
                        CNrt$myvar,LMA$myvar)

names(all_predictors) <- c("lon","lat","z","Tg","PPFD","vpd",
                           "alpha","fAPAR","age","CNrt","LMA")

Tg_df <- all_predictors[,c("lon","lat","z","Tg")]
PPFD_df <- all_predictors[,c("lon","lat","z","PPFD")]
vpd_df <- all_predictors[,c("lon","lat","z","vpd")]
alpha_df <- all_predictors[,c("lon","lat","z","alpha")]
fAPAR_df <- all_predictors[,c("lon","lat","z","fAPAR")]
age_df <- all_predictors[,c("lon","lat","z","age")]
CNrt_df <- all_predictors[,c("lon","lat","z","CNrt")]
LMA_df <- all_predictors[,c("lon","lat","z","LMA")]

#now, apply gwr to extract site predictors' value
NPP_Forest$Tg <- NA
NPP_Forest$PPFD <- NA
NPP_Forest$vpd <- NA
NPP_Forest$alpha <- NA
NPP_Forest$fAPAR <- NA
NPP_Forest$age <- NA
NPP_Forest$CNrt <- NA
NPP_Forest$LMA <- NA
NPP_Forest$z <- NPP_Forest$elv
a <- 1.5 # which degree (distance) of grid when interpolating gwr from global grids
#Extract Tg, PPFD, vpd, alpha,fAPAR,age,CNrt,LMA, max-vcmax25
for (i in 1:nrow(NPP_Forest)) {
  tryCatch({
    #Tg
    Tg_global <- na.omit(Tg_df)
    NRE_part <- subset(Tg_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("Tg")] <- (gwr(Tg ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #ppfd
    PPFD_global <- na.omit(PPFD_df)
    NRE_part <- subset(PPFD_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("PPFD")] <- (gwr(PPFD ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #vpd
    vpd_global <- na.omit(vpd_df)
    NRE_part <- subset(vpd_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("vpd")] <- (gwr(vpd ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #alpha
    alpha_global <- na.omit(alpha_df)
    NRE_part <- subset(alpha_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("alpha")]  <- (gwr(alpha ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #fAPAR
    fAPAR_global <- na.omit(fAPAR_df)
    NRE_part <- subset(fAPAR_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("fAPAR")]<- (gwr(fAPAR ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #age
    age_global <- na.omit(age_df)
    NRE_part <- subset(age_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("age")]  <- (gwr(age ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #CNrt
    CNrt_global <- na.omit(CNrt_df)
    NRE_part <- subset(CNrt_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("CNrt")]  <- (gwr(CNrt ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #LMA
    LMA_global <- na.omit(LMA_df)
    NRE_part <- subset(LMA_global,lon>(NPP_Forest[i,"lon"]-a)&lon<(NPP_Forest[i,"lon"]+a)&
                         lat>(NPP_Forest[i,"lat"]-a)&lat<(NPP_Forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_Forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_Forest[i,c("LMA")]  <- (gwr(LMA ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
  }, error=function(e){})} 

summary(NPP_Forest)
NPP_Forest$vpd[NPP_Forest$vpd<0] <- NA
NPP_Forest$age[NPP_Forest$age<0] <- NA

#load statistical model coefficients
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

#now, using several statistical models to predict npp, anpp, npp.leaf....
#now, using several statistical models to predict npp, anpp, npp.leaf....
NPP_Forest$pred_npp <- NPP_Forest$pred_gpp_c3 * (1/(1 + exp(-(summary(mod_tnpp)$coefficients[1,1] + 
                                                                summary(mod_tnpp)$coefficients[2,1] * log(NPP_Forest$CNrt) +
                                                                summary(mod_tnpp)$coefficients[3,1] * log(NPP_Forest$age) + 
                                                                summary(mod_tnpp)$coefficients[4,1] * NPP_Forest$fAPAR))))

NPP_Forest$pred_anpp <- NPP_Forest$pred_gpp_c3 * (1/(1 + exp(-(summary(mod_anpp)$coefficients[1,1]+
                                                                 summary(mod_anpp)$coefficients[2,1] * log(NPP_Forest$CNrt) +
                                                                 summary(mod_anpp)$coefficients[3,1] * log(NPP_Forest$age) + 
                                                                 summary(mod_anpp)$coefficients[4,1] * NPP_Forest$fAPAR))))

NPP_Forest$pred_bnpp <- NPP_Forest$pred_npp - NPP_Forest$pred_anpp

NPP_Forest$pred_lnpp <- NPP_Forest$pred_anpp * (1/(1 + exp(-(summary(mod_lnpp)$coefficients[1,1]+
                                                               summary(mod_lnpp)$coefficients[2,1]* log(NPP_Forest$PPFD) +
                                                               summary(mod_lnpp)$coefficients[3,1] * (NPP_Forest$Tg) + 
                                                               summary(mod_lnpp)$coefficients[4,1] * log(NPP_Forest$vpd)))))

NPP_Forest$pred_wnpp <- NPP_Forest$pred_anpp - NPP_Forest$pred_lnpp

#use rsofun - site-species
#NPP_Forest$pred_leafnc <- (0.0162/0.5) + (0.0039/0.5) * NPP_Forest$max_vcmax25/NPP_Forest$LMA
NPP_Forest$pred_leafnc <- (summary(n1)$coefficients[1,1]/0.46) + (summary(n1)$coefficients[2,1]/0.46) * NPP_Forest$max_vcmax25/NPP_Forest$LMA

NPP_Forest$pred_lnf <- NPP_Forest$pred_lnpp*NPP_Forest$pred_leafnc

#summary(NPP_Forest$CN_wood_final)
NPP_Forest$pred_wnf <- NPP_Forest$pred_wnpp/100
NPP_Forest$pred_bnf <- NPP_Forest$pred_bnpp/94
#root mean was obtained from median of NPP forest

NPP_Forest$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] + summary(nre_model)$coefficients[2,1] *NPP_Forest$Tg + summary(nre_model)$coefficients[3,1] * log(NPP_Forest$vpd)))))

hist(NPP_Forest$pred_nre)

NPP_Forest$pred_nuptake <- NPP_Forest$pred_lnf*(1-NPP_Forest$pred_nre) + NPP_Forest$pred_wnf +NPP_Forest$pred_bnf  

NPP_Forest_gcme <- subset(NPP_Forest,method=="Total Nuptake reported in GCME, at ambient condition")
dim(NPP_Forest_gcme)
My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 20))

#read org data again, to get info of fertilizations
df <- read.csv("/Users/yunpeng/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv")
gcme_nup <- subset(df,Data_type=="N_uptake")
dim(gcme_nup)
#firstly, for co2
gcme_nup_co2 <- gcme_nup[,c("co2_a","co2_e","exp_nam","treatment")]
gcme_nup_co2$co2_a <- as.numeric(gcme_nup_co2$co2_a)
gcme_nup_co2$co2_e <- as.numeric(gcme_nup_co2$co2_e)
summary(gcme_nup_co2)
gcme_nup_co2 <- aggregate(gcme_nup_co2,by=list(gcme_nup_co2$exp_nam,gcme_nup_co2$treatment), FUN=mean, na.rm=TRUE) #site-mean

gcme_nup_co2_final <- gcme_nup_co2[,c(1:4)]
names(gcme_nup_co2_final) <- c("exp_nam","treatment","co2_a","co2_e")
#merge

NPP_Forest_gcme_co2 <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","treatment"),all.x=TRUE), 
                   list(NPP_Forest_gcme,gcme_nup_co2_final))
dim(NPP_Forest_gcme_co2)

#now, n addition --> have a look
gcme_nup_n <- gcme_nup[,c("nfertQ_a","nfertQ_e","exp_nam","treatment")]
gcme_nup_n <- subset(gcme_nup_n,treatment!="c" & treatment!="w"&treatment!="cw")
gcme_nup_n
#michgan data looks strange --> but it seems that ambient = 0, elevated = 60
#gcme_nup_n$nfertQ_a[gcme_nup_n$exp_nam=="Michigan_UNDERC_bog"] <- 0
#gcme_nup_n$nfertQ_a[gcme_nup_n$exp_nam=="Michigan_UNDERC_intermFen"]<- 0
#gcme_nup_n$nfertQ_a[gcme_nup_n$exp_nam=="Michigan_UNDERC_richFen"]<- 0
#gcme_nup_n$nfertQ_e[gcme_nup_n$exp_nam=="Michigan_UNDERC_bog"] <- 60
#gcme_nup_n$nfertQ_e[gcme_nup_n$exp_nam=="Michigan_UNDERC_intermFen"]<- 60
#gcme_nup_n$nfertQ_e[gcme_nup_n$exp_nam=="Michigan_UNDERC_richFen"]<- 60

gcme_nup_n <- aggregate(gcme_nup_n,by=list(gcme_nup_n$exp_nam,gcme_nup_n$treatment), FUN=mean, na.rm=TRUE) #site-mean

gcme_nup_n_final <- gcme_nup_n[,c(1:4)]
names(gcme_nup_n_final) <- c("exp_nam","treatment","nitrogen_a","nitrogen_e")
#merge

NPP_Forest_gcme_co2_n <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","treatment"),all.x=TRUE), 
                             list(NPP_Forest_gcme_co2,gcme_nup_n_final))
head(NPP_Forest_gcme_co2_n)

hist(NPP_Forest_gcme_co2_n$obs_nuptake)
hist(NPP_Forest_gcme_co2_n$elevated)

NPP_Forest_gcme_co2_n$N_GPP_a <- NPP_Forest_gcme_co2_n$obs_nuptake/NPP_Forest_gcme_co2_n$pred_gpp_c3
NPP_Forest_gcme_co2_n$N_GPP_e <- NPP_Forest_gcme_co2_n$elevated/NPP_Forest_gcme_co2_n$pred_gpp_c3

NPP_Forest_gcme_co2_n$nuptake_change <- (NPP_Forest_gcme_co2_n$elevated-NPP_Forest_gcme_co2_n$obs_nuptake)/NPP_Forest_gcme_co2_n$obs_nuptake
NPP_Forest_gcme_co2_n$N_GPP_change <- (NPP_Forest_gcme_co2_n$N_GPP_e-NPP_Forest_gcme_co2_n$N_GPP_a)/NPP_Forest_gcme_co2_n$N_GPP_a

NPP_Forest_gcme_co2_n$treatment[NPP_Forest_gcme_co2_n$treatment=="f2"] <- "f"

NPP_Forest_gcme_co2_n %>% group_by(exp_nam) %>% 
  summarise (number=n())

#why it is different? let's using treatment firstly?
NPP_Forest_gcme_co2_n %>% group_by(treatment) %>% summarise (number=n())
NPP_Forest_gcme_co2_n %>% group_by(factors) %>%summarise (number=n())

boxplot(nuptake_change~treatment,data=NPP_Forest_gcme_co2_n,
        main="Change of total N uptake (%)",
        xlab="treatment", ylab="N uptake change (%)",ylim = c(0, 20),
        names=c("CO2","CO2 + N addition","CO2 + warming",
                "N addition","N addition + warming","warming"))

boxplot(N_GPP_change~treatment,data=NPP_Forest_gcme_co2_n,
        main="Change of N demand of GPP (%)",
        xlab="treatment", ylab="N demand of GPP change (%)",ylim = c(0, 20),
        names=c("CO2","CO2 + N addition","CO2 + warming",
                "N","N addition + warming","warming"))

#merged with measured npp
df <- read.csv("/Users/yunpeng/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv")
gcme_npp <- subset(df,Data_type=="NPP")
dim(gcme_npp)
#firstly, for co2
gcme_npp <- gcme_npp[,c("ambient","elevated","exp_nam","treatment")]
gcme_npp$elevated <- as.numeric(gcme_npp$elevated)
gcme_npp$ambient <- as.numeric(gcme_npp$ambient)
gcme_npp <- aggregate(gcme_npp,by=list(gcme_npp$exp_nam,gcme_npp$treatment), FUN=mean, na.rm=TRUE) #site-mean

gcme_npp_final <- gcme_npp[,c(1:4)]
names(gcme_npp_final) <- c("exp_nam","treatment","ambient_npp","elevated_npp")
#merge

NPP_Forest_gcme_co2_n_npp <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","treatment"),all.x=TRUE), 
                             list(NPP_Forest_gcme_co2_n,gcme_npp_final))
dim(NPP_Forest_gcme_co2_n_npp)
summary(NPP_Forest_gcme_co2_n_npp)

NPP_Forest_gcme_co2_n_npp$N_NPP_a <- NPP_Forest_gcme_co2_n_npp$obs_nuptake/NPP_Forest_gcme_co2_n_npp$ambient_npp
NPP_Forest_gcme_co2_n_npp$N_NPP_e <- NPP_Forest_gcme_co2_n_npp$elevated/NPP_Forest_gcme_co2_n_npp$elevated_npp

NPP_Forest_gcme_co2_n_npp$N_NPP_change <- (NPP_Forest_gcme_co2_n_npp$N_NPP_e-NPP_Forest_gcme_co2_n_npp$N_NPP_a)/NPP_Forest_gcme_co2_n_npp$N_NPP_a

boxplot(N_NPP_change~treatment,data=NPP_Forest_gcme_co2_n_npp,
        main="Change of N demand of NPP (%)",
        xlab="treatment", ylab=" ",
        names=c("co2","N addition","N addition+warming","warming"))
dim(subset(NPP_Forest_gcme_co2_n_npp,is.na(N_NPP_change)==FALSE))


#merged with measured gpp
df <- read.csv("/Users/yunpeng/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv")
gcme_gpp <- subset(df,Data_type=="GPP")
dim(gcme_gpp)
#firstly, for co2
gcme_gpp <- gcme_gpp[,c("ambient","elevated","exp_nam","treatment")]
gcme_gpp$elevated <- as.numeric(gcme_gpp$elevated)
gcme_gpp$ambient <- as.numeric(gcme_gpp$ambient)
gcme_gpp <- aggregate(gcme_gpp,by=list(gcme_gpp$exp_nam,gcme_gpp$treatment), FUN=mean, na.rm=TRUE) #site-mean

gcme_gpp_final <- gcme_gpp[,c(1:4)]
names(gcme_gpp_final) <- c("exp_nam","treatment","ambient_gpp","elevated_gpp")
#merge

NPP_Forest_gcme_co2_n_npp_gpp <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","treatment"),all.x=TRUE), 
                                   list(NPP_Forest_gcme_co2_n_npp,gcme_gpp_final))
dim(NPP_Forest_gcme_co2_n_npp_gpp)
summary(NPP_Forest_gcme_co2_n_npp_gpp)
