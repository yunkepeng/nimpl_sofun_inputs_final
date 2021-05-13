rm(list=ls())
load(file = "/Users/yunpeng/data/NPP_final/Forest_site_simulation.Rdata")

devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/rsofun/")
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

#input the current complete version of NPP_forest
#before start everything - convert the measurement year before 1980 to 1980-1989 (which is consistent with what we set in climate forcing), so that we can run them sucessfully in rsofun later on.

#For how to collect them. see L1-714 of "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_site_orig.R" - sitename and sitename_fapar already been well checked
NPP_Forest <- read.csv("/Users/yunpeng/data/NPP_final/NPP_Forest.csv")
NPP_Forest$year_start <-NPP_Forest$Begin_year
NPP_Forest$year_end <-NPP_Forest$End_year

NPP_Forest$year_start[NPP_Forest$Begin_year<=1980] <- 1980
NPP_Forest$year_end[NPP_Forest$End_year<=1980] <- 1989

####now, input forcing data from two times simulation
forcing_df <- list.files("/Users/yunpeng/data/NPP_final/reprocessing_climates/",full.names = T)
length(forcing_df)

fapar_df <- list.files("/Users/yunpeng/data/NPP_final/reprocessing_fpar/",full.names = T)
length(fapar_df)

#1. fapar - input

#1. fpar - check missing data - and also, input all years fapar (2001-2015), which will be selected in measurement year only later on 
for (i in 1:length(fapar_df)){
  df1 <- read.csv(fapar_df[i])
  df1$date <- as.Date(df1$date)
  df1 <- df1[!(format(df1$date,"%m") == "02" & format(df1$date, "%d") == "29"), , drop = FALSE]
  df2 <- df1[,c("date","modisvar_filled")]
  assign(substr(sub('.*daily_', '', fapar_df[i]),1,nchar(sub('.*daily_', '', fapar_df[i]))-4), df2) 
}

#check fapar missing data
for (i in 1:nrow(NPP_Forest)){
  NPP_Forest$forcing_avil[i] <- exists(paste(NPP_Forest$sitename_fpar[i]))
}

na_fapar <- (subset(NPP_Forest,forcing_avil=="FALSE" ))
dim(na_fapar)
na_fapar$sitename

library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(na_fapar$lon,na_fapar$lat, col="red", pch=16,cex=1)
#3 samples were missing fapar in the edge, which was expected.

#check climate forcing missing data
empty_vec <- c()

#check existed climate files
for (i in 1:(length(forcing_df))){
  empty_vec[i] <- as.numeric(gsub("[^0-9]", "",  forcing_df[i]))
}

diff <- setdiff(1:935, empty_vec)
diff
NPP_Forest[diff,]

points(NPP_Forest$lon[338],NPP_Forest$lat[338], col="blue", pch=16,cex=1)
points(NPP_Forest$lon[390],NPP_Forest$lat[390], col="blue", pch=16,cex=1)

#totally 5 sites were missing:
na_fapar$sitename
NPP_Forest[diff,]$sitename
all_na_points <- c(na_fapar$sitename,NPP_Forest[diff,]$sitename)
#NPP_F556, NPP_F697, NPP_F700 (due to fapar) and NPP_F338, NPP_F390 (due to climate forcing)



#2. forcing - combing fapar and climates into a df.

for (i in 1:length(forcing_df)){
  tryCatch({
  df1 <- read.csv(forcing_df[i])
  df1$date <- as.Date(df1$date)
  
  sitename_climate <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$sitename
  sitename_fapar <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$sitename_fpar
  # 5 points were missing - let's clarify them firstly
  if (sitename_climate %in% all_na_points){
    print (sitename_climate)
    print ("this site is not available")
  } else {
    fapar <- (eval(parse(text=sitename_fapar)))
    fapar$Year <- year(fapar$date)
    
    yr_start <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$year_start
    yr_end <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$year_end
    
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
}, error=function(e){})} 



#3. rsofun to predict gpp
df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286)

NPP_Forest$whc = 170

NPP_Forest$pred_gpp_c3 <- NA
#NPP_Forest$pred_gpp_c4 <- NA
NPP_Forest$max_vcmax25_c3 <- NA
#NPP_Forest$max_vcmax25_c4 <- NA


#using rsofun
for (i in 1:nrow(NPP_Forest)) {
  tryCatch({
    #c3
    forcing <- (eval(parse(text=(paste("final",NPP_Forest$sitename[i],sep="_")))))
    modlist <- run_pmodel_f_bysite( 
      NPP_Forest$sitename[i], 
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
        firstyeartrend = NPP_Forest$year_start[i],
        nyeartrend = NPP_Forest$year_end[i]-NPP_Forest$year_start[i]+1), 
      siteinfo = NPP_Forest[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    
    pred_gpp_list <- modlist %>% mutate(ymonth = month(date),yday = day(date)) %>% group_by(ymonth, yday) %>% summarise(gpp = mean(gpp, na.rm = TRUE))
    max_vcmax25 <- max(modlist$vcmax25)*1000000
    
    NPP_Forest[i,c("pred_gpp_c3")] <- sum(pred_gpp_list$gpp)
    NPP_Forest[i,c("max_vcmax25_c3")] <- max_vcmax25
  }, error=function(e){})} 

#this sites have no gpp data - must because their fapar in n_focal = 0, we need to fill them by alternatively applying n_focal = 1, then 2...
subset(NPP_Forest,pred_gpp_c3=="NaN")$sitename_fpar

#we have newly interpolate their fapar primarily based on n_focal = 1, then n_focal = 2, and saved it in "/Users/yunpeng/data/forest_npp/reprocessing_fpar_raw/"
#the code of this is available at L90-110 in forest/Reprocessing_fpar_climates_forest.R

#now, reprocessing such values - by updating such fapar 
fapar_df_new <- list.files("~/data/forest_npp/reprocessing_fpar_raw/",full.names = T)

for (i in 1:(length(fapar_df_new)-1)){
  df1 <- read.csv(fapar_df_new[i])
  df1$date <- as.Date(df1$date)
  df1 <- df1[!(format(df1$date,"%m") == "02" & format(df1$date, "%d") == "29"), , drop = FALSE]
  df2 <- df1[,c("date","modisvar_filled")]
  assign(substr(sub('.*daily_', '', fapar_df_new[i]),1,nchar(sub('.*daily_', '', fapar_df_new[i]))-4), df2) 
}


#2. forcing - combing fapar and climates into a dataframe.
for (i in 1:length(forcing_df)){
  df1 <- read.csv(forcing_df[i])
  df1$date <- as.Date(df1$date)
  
  sitename_climate <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$sitename
  sitename_fapar <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$sitename_fpar
  # 5 points were missing - let's clarify them firstly
  if (sitename_climate %in% all_na_points){
    print (sitename_climate)
    print ("this site is not available")
  } else {
    fapar <- (eval(parse(text=sitename_fapar)))
    fapar$Year <- year(fapar$date)
    
    yr_start <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$year_start
    yr_end <- subset(NPP_Forest,NPP_Forest$sitename == df1$sitename[1])$year_end
    
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
}

#3. rsofun to predict gpp
df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286)

NPP_Forest$whc = 170

NPP_Forest$pred_gpp_c3 <- NA
#NPP_Forest$pred_gpp_c4 <- NA
NPP_Forest$max_vcmax25_c3 <- NA
#NPP_Forest$max_vcmax25_c4 <- NA


#using rsofun
for (i in 1:nrow(NPP_Forest)) {
  tryCatch({
    #c3
    forcing <- (eval(parse(text=(paste("final",NPP_Forest$sitename[i],sep="_")))))
    modlist <- run_pmodel_f_bysite( 
      NPP_Forest$sitename[i], 
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
        firstyeartrend = NPP_Forest$year_start[i],
        nyeartrend = NPP_Forest$year_end[i]-NPP_Forest$year_start[i]+1), 
      siteinfo = NPP_Forest[i,], 
      forcing, 
      df_soiltexture, 
      params_modl = params_modl, 
      makecheck = TRUE)
    
    pred_gpp_list <- modlist %>% mutate(ymonth = month(date),yday = day(date)) %>% group_by(ymonth, yday) %>% summarise(gpp = mean(gpp, na.rm = TRUE))
    max_vcmax25 <- max(modlist$vcmax25)*1000000
    
    NPP_Forest[i,c("pred_gpp_c3")] <- sum(pred_gpp_list$gpp)
    NPP_Forest[i,c("max_vcmax25_c3")] <- max_vcmax25
  }, error=function(e){})} 

subset(NPP_Forest,is.na(pred_gpp_c3)==TRUE)$sitename_fpar

plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(subset(NPP_Forest,is.na(pred_gpp_c3)==TRUE)$lon,subset(NPP_Forest,is.na(pred_gpp_c3)==TRUE)$lat, col="red", pch=16,cex=1)
#these points were missing, either due to fapar or climate forcing missing

NPP_Forest_all_flux <- NPP_Forest
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
#plot lma for a second
coordinates(LMA) <- ~lon+lat 
gridded(LMA) <- TRUE

r4 <- raster(LMA, "myvar") 
plot(r4)


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

newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(subset(NPP_Forest,is.na(LMA)==TRUE)$lon,subset(NPP_Forest,is.na(LMA)==TRUE)$lat, col="red", pch=16,cex=1)
#these sites have no data for LMA, CNrt age and fAPAR.

#now, using several statistical models to predict npp, anpp, npp.leaf....
NPP_Forest$pred_npp <- NPP_Forest$pred_gpp_c3 * (1/(1 + exp(-(-0.36075 * log(NPP_Forest$CNrt) +
                                                                -0.16213 * log(NPP_Forest$age) + 
                                                                0.72793 * NPP_Forest$fAPAR+ 0.57014))))

NPP_Forest$pred_anpp <- NPP_Forest$pred_gpp_c3 * (1/(1 + exp(-(-0.55151 * log(NPP_Forest$CNrt) +
                                                              -0.20050 * log(NPP_Forest$age) + 
                                                              1.06611 * NPP_Forest$fAPAR+ 0.35817))))

NPP_Forest$pred_bnpp <- NPP_Forest$pred_npp - NPP_Forest$pred_anpp

NPP_Forest$pred_lnpp <- NPP_Forest$pred_anpp * (1/(1 + exp(-(0.97093* log(NPP_Forest$PPFD) +
                                                               0.06453 * (NPP_Forest$Tg) + 
                                                               -0.80397 * log(NPP_Forest$vpd) + -7.47165))))

NPP_Forest$pred_wnpp <- NPP_Forest$pred_anpp - NPP_Forest$pred_lnpp

#use rsofun - site-species
#NPP_Forest$pred_leafnc <- (0.0162/0.5) + (0.0039/0.5) * NPP_Forest$max_vcmax25/NPP_Forest$LMA
NPP_Forest$pred_leafnc <- (0.01599/0.46) + (0.005992/0.46) * NPP_Forest$max_vcmax25/NPP_Forest$LMA

NPP_Forest$pred_lnf <- NPP_Forest$pred_lnpp*NPP_Forest$pred_leafnc

#summary(NPP_Forest$CN_wood_final)
NPP_Forest$pred_wnf <- NPP_Forest$pred_wnpp/100

hist(NPP_Forest$CN_root_final)
summary(NPP_Forest$CN_root_final)

NPP_Forest$pred_bnf <- NPP_Forest$pred_bnpp/122

NPP_Forest$CN_leaf_final[NPP_Forest$CN_leaf_final>100] <- NA
NPP_Forest$CN_leaf_org[NPP_Forest$CN_leaf_org>100] <- NA

#the difference between cn_leaf_org and cn_leaf_final is that the previous one has not include "repeated merge" (with 36 less sites)

#correct new dataset's rep_info
NPP_Forest$rep_info[is.na(NPP_Forest$rep_info)==TRUE] <- ""

#median of wood cn and root cn
#using median of root = 94
#using median of wood =100
# see below
summary(read.csv("/Users/yunpeng/data/CN_wood/wood_cn.csv")$OrigValueStr)
#summary(NPP_Forest2$CN_root_final)
#Remove rep and rep2 (rep1 and rep2 are paired repetation so we can safely remover rep2 since it is ForC, while rep1 is from Sara Vicca)
NPP_Forest2 <- subset(NPP_Forest,rep_info!="rep" & rep_info!="rep2" &file!="NPP_Schulze")
#Remove Schulze - their c/n unit is kg/ha - which is quite strange - we don't know how many species they considered for leaf C/N! See hist of their lnf_obs_org. Also, the reference is too old, which is not reliable when just citing them in a book
summary(NPP_Forest2$CN_root_final)

NPP_Forest2_sitemean <- aggregate(NPP_Forest2,by=list(NPP_Forest2$lon,NPP_Forest2$lat,NPP_Forest2$z), FUN=mean, na.rm=TRUE) #site-mean

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 20))

#check
#analyse_modobs2(forest_site2,"pred_gpp", "GPP",type = "points")
ggplot(data=NPP_Forest2, aes(x=pred_gpp_c3, y=GPP)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(GPP~pred_gpp_c3,NPP_Forest2))


#analyse_modobs2(forest_site2,"pred_npp", "TNPP_1",type = "points")
ggplot(data=NPP_Forest2, aes(x=pred_npp, y=TNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(TNPP_1~pred_npp,NPP_Forest2))

#analyse_modobs2(forest_site2,"pred_anpp", "ANPP_2",type = "points")
ggplot(data=NPP_Forest2, aes(x=pred_anpp, y=ANPP_2)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(ANPP_2~pred_anpp,NPP_Forest2))

#analyse_modobs2(forest_site2,"pred_lnpp", "NPP.foliage",type = "points")
ggplot(data=NPP_Forest2, aes(x=pred_lnpp, y=NPP.foliage)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(NPP.foliage~pred_lnpp,NPP_Forest2))

#analyse_modobs2(forest_site,"pred_wnpp", "NPP.wood",type = "points")
ggplot(data=NPP_Forest2, aes(x=pred_wnpp, y=NPP.wood)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(NPP.wood~pred_wnpp,NPP_Forest2))

#analyse_modobs2(forest_site2,"pred_bnpp", "BNPP_1",type = "points")
ggplot(data=NPP_Forest2, aes(x=pred_bnpp, y=BNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(BNPP_1~pred_bnpp,NPP_Forest2))

#analyse_modobs2(forest_site,"pred_lnf", "lnf_obs",type = "points") 
ggplot(data=NPP_Forest2, aes(x=pred_lnf, y=lnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(lnf_obs_final~pred_lnf,NPP_Forest2))

#without second interpolation --> use this!
ggplot(data=NPP_Forest2, aes(x=pred_lnf, y=lnf_obs_org)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(lnf_obs_org~pred_lnf,NPP_Forest2))

#wnf - assuming constant wood/cn = 97

ggplot(data=NPP_Forest2, aes(x=pred_wnf, y=wnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(wnf_obs_final~pred_wnf,NPP_Forest2))

#bnf - assuming constant root/cn = 122
ggplot(data=NPP_Forest2, aes(x=pred_bnf, y=bnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(bnf_obs_final~pred_bnf,NPP_Forest2))

NPP_Forest2 %>% group_by(file) %>% summarise(number = n())

#leaf cn
#(9) leafcn
#check leaf c/n
SP_input <- read.csv(file="/Users/yunpeng/data/leaf_traits/combined_leaf_traits.csv") #new one 
SP_input2 <- SP_input[,c("lat","lon","z","Vcmax25","narea","lma")]
sitemean <- aggregate(SP_input2,by=list(SP_input2$lon,SP_input2$lat), FUN=mean, na.rm=TRUE) 
dim(sitemean)

sitemean$pred_leafn <- (0.01599) + (0.005992)* sitemean$Vcmax25/sitemean$lma
sitemean$obs_leafn <- sitemean$narea/sitemean$lma

ggplot(data=sitemean, aes(x=pred_leafn, y=obs_leafn)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(obs_leafn~pred_leafn,sitemean))

#nuptake
NPP_Forest2$pred_nuptake <- NPP_Forest2$pred_lnf + NPP_Forest2$pred_bnf + NPP_Forest2$pred_wnf
NPP_Forest2$obs_nuptake <- NPP_Forest2$lnf_obs_org + NPP_Forest2$bnf_obs_final + NPP_Forest2$wnf_obs_final

ggplot(data=NPP_Forest2, aes(x=pred_nuptake, y=obs_nuptake)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(pred_nuptake~obs_nuptake,NPP_Forest2))

#nre
#check mean of NRE
NRE_Du <- read.csv(file="~/data/NRE_various/NRE_Du/NRE_Du.csv")
NRE_Dong <- read.csv(file="~/data/NRE_various/NRE_Deng/NRE_Deng.csv")

NRE_Du_df <- NRE_Du[,c("lon","lat","NRE","MAT","MAP")]
NRE_Du_df <- aggregate(NRE_Du_df,by=list(NRE_Du_df$lon,NRE_Du_df$lat), FUN=mean, na.rm=TRUE) #site-mean
NRE_Du_df <- NRE_Du_df[,c(3:7)]
head(NRE_Du_df)
dim(NRE_Du_df)

NRE_Dong_df <- NRE_Dong[,c("Longitude","Latitude","NRE.nitrogen.resorption.efficiency.","MAT","MAP")]
names(NRE_Dong_df) <- c("lon","lat","NRE","MAT","MAP")
head(NRE_Dong_df)
NRE_Dong_df <- aggregate(NRE_Dong_df,by=list(NRE_Dong_df$lon,NRE_Dong_df$lat), FUN=mean, na.rm=TRUE) #site-mean
NRE_Dong_df <- NRE_Dong_df[,c(3:7)]
dim(NRE_Dong_df)


NRE_Dong_df$source <- "Dong"
NRE_Du_df$source <- "Du"
NRE_df <- rbind(NRE_Du_df,NRE_Dong_df)
summary(NRE_df)

#check repeated data, and remove 6 repeated points from Du et al. paper
NRE_df$repeated <- duplicated(NRE_df[,c("lon","lat")])
summary(NRE_df$repeated)
NRE_df <- subset(NRE_df,repeated==FALSE)

#project data
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)

points(NRE_df$lon,NRE_df$lat, col="red", pch=16,cex=1)

#3. add elevation in this df, based on ingtestr 
siteinfo <- NRE_df[,c("lon","lat")] # present x and y separately
siteinfo$date_start <- lubridate::ymd(paste0(1982, "-01-01"))
siteinfo$date_end <- lubridate::ymd(paste0(2011, "-12-31"))
siteinfo$sitename <- paste0("s", 1:nrow(siteinfo),sep="")
siteinfo <- as_tibble(siteinfo)

devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/ingestr/")


df_etopo <- ingest(
  siteinfo,
  source = "etopo1",
  dir = "~/data/etopo/" 
)

NRE_df$elevation <- as.numeric(as.data.frame(df_etopo$data))
subset(NRE_df,elevation<0)
#Some grids > 0, lets' assume -3062 as NA, and others as 0 firstly?
NRE_df$elevation[NRE_df$elevation< -50] <- NA
NRE_df$elevation[NRE_df$elevation< 0] <- 0

summary(NRE_df)

#for nre
names(NRE_df) <- c("lon","lat","NRE","MAT","MAP","source","repeated","z")
NRE_df$Tg <- NA
NRE_df$vpd <- NA
a <- 1.5

for (i in 1:nrow(NRE_df)) {
  tryCatch({
    #Tg
    Tg_global <- na.omit(Tg_df)
    NRE_part <- subset(Tg_global,lon>(NRE_df[i,1]-a)&lon<(NRE_df[i,1]+a)&
                         lat>(NRE_df[i,2]-a)&lat<(NRE_df[i,2]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NRE_df[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NRE_df[i,c("Tg")] <- (gwr(Tg ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #vpd
    vpd_global <- na.omit(vpd_df)
    NRE_part <- subset(vpd_global,lon>(NRE_df[i,1]-a)&lon<(NRE_df[i,1]+a)&
                         lat>(NRE_df[i,2]-a)&lat<(NRE_df[i,2]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NRE_df[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NRE_df[i,c("vpd")] <- (gwr(vpd ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
  }, error=function(e){})} 

NRE_df$pred_nre <- NA
NRE_df$vpd[NRE_df$vpd<0] <- NA
NRE_df$pred_nre <- (1/(1+exp(-(-0.064460 *NRE_df$Tg + 0.402850 * log(NRE_df$vpd) + 1.368935))))

NRE_df$NRE <- NRE_df$NRE/100

ggplot(data=NRE_df, aes(x=pred_nre, y=NRE)) + xlim(c(0.25,1))+ylim(c(0.25,1))+
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme

summary(lm(NRE~pred_nre,NRE_df))

#now, newly adding fluxnet gpp sites
fluxnet_site <- ingestr::siteinfo_fluxnet2015
fluxnet_site$lat <- as.numeric(format(round(fluxnet_site$lat, 2), nsmall = 2))
fluxnet_site$lon <- as.numeric(format(round(fluxnet_site$lon, 2), nsmall = 2))

GMD_GPP <- read.csv(file="~/data/gpp_gmd/gpp_pmodel_fluxnet2015_stocker19gmd_spatial.csv")
GMD_GPP <- subset(GMD_GPP,setup=="FULL")

fluxnet_GPP <- merge(fluxnet_site,GMD_GPP,by=c("sitename"),all.x=TRUE)

fluxnet_GPP2 <-fluxnet_GPP[,c("lon","lat","sitename","gpp")]
fluxnet_GPP2 <- na.omit(fluxnet_GPP2)

NPP_site <- NPP_Forest[,c("lon","lat","z","site")]
NPP_site$lat <- as.numeric(format(round(NPP_site$lat, 2), nsmall = 2))
NPP_site$lon <- as.numeric(format(round(NPP_site$lon, 2), nsmall = 2))


NPP_test <- merge(NPP_site,fluxnet_GPP2,by=c("lon","lat"),all.x=TRUE)
NPP_test<- na.omit(NPP_test)

fluxnet_gpp_final <- aggregate(NPP_test,by=list(NPP_test$site), mean,na.rm=TRUE)
fluxnet_gpp_final <- fluxnet_gpp_final[,c(1,7)]
names(fluxnet_gpp_final) <- c("site","flux_gpp")
fluxnet_gpp_final

#now, interpolating 2nd version (with the inclusion of fluxnet predicted gpp)
NPP_Forest_by_flux <- merge(NPP_Forest,fluxnet_gpp_final,by=c("site"),all.x=TRUE)

NPP_Forest_by_flux$new_gpp <- NA

for (i in 1:nrow(NPP_Forest_by_flux)){
  if (is.na(NPP_Forest_by_flux$flux_gpp[i]) == TRUE){ 
    NPP_Forest_by_flux$new_gpp[i] <- NPP_Forest_by_flux$pred_gpp_c3[i]
  } else {
    NPP_Forest_by_flux$new_gpp[i] <- NPP_Forest_by_flux$flux_gpp[i]
  }
}
summary(NPP_Forest_by_flux$flux_gpp)

#now, repeat


#now, using several statistical models to predict npp, anpp, npp.leaf....

NPP_Forest_by_flux$pred_npp <- NPP_Forest_by_flux$new_gpp * (1/(1 + exp(-(-0.36075 * log(NPP_Forest_by_flux$CNrt) +
                                                                            -0.16213 * log(NPP_Forest_by_flux$age) + 
                                                                            0.72793 * NPP_Forest_by_flux$fAPAR + 0.57014))))

NPP_Forest_by_flux$pred_anpp <- NPP_Forest_by_flux$new_gpp * (1/(1 + exp(-(-0.55151 * log(NPP_Forest_by_flux$CNrt) +
                                                              -0.20050 * log(NPP_Forest_by_flux$age) + 
                                                                1.06611 * NPP_Forest_by_flux$fAPAR+ 0.35817))))


NPP_Forest_by_flux$pred_bnpp <- NPP_Forest_by_flux$pred_npp - NPP_Forest_by_flux$pred_anpp

NPP_Forest_by_flux$pred_lnpp <- NPP_Forest_by_flux$pred_anpp * (1/(1 + exp(-(0.97093* log(NPP_Forest_by_flux$PPFD) +
                                                                               0.06453 * (NPP_Forest_by_flux$Tg) + 
                                                               -0.80397 * log(NPP_Forest_by_flux$vpd) + -7.47165))))

NPP_Forest_by_flux$pred_wnpp <- NPP_Forest_by_flux$pred_anpp - NPP_Forest_by_flux$pred_lnpp

#use rsofun - site-species
#NPP_Forest_by_flux$pred_leafnc <- (0.0162/0.5) + (0.0039/0.5) * NPP_Forest_by_flux$max_vcmax25/NPP_Forest_by_flux$LMA
NPP_Forest_by_flux$pred_leafnc <- (0.01599/0.46) + (0.005992/0.46) * NPP_Forest_by_flux$max_vcmax25/NPP_Forest_by_flux$LMA


NPP_Forest_by_flux$pred_lnf <- NPP_Forest_by_flux$pred_lnpp*NPP_Forest_by_flux$pred_leafnc

NPP_Forest_by_flux$pred_wnf <- NPP_Forest_by_flux$pred_wnpp/97

NPP_Forest_by_flux$pred_bnf <- NPP_Forest_by_flux$pred_bnpp/122

#correct new dataset's rep_info
NPP_Forest_by_flux2 <- subset(NPP_Forest_by_flux,rep_info!="rep" & rep_info!="rep1"& file!="NPP_Schulze")

NPP_Forest_by_flux2_sitemean <- aggregate(NPP_Forest_by_flux2,by=list(NPP_Forest_by_flux2$lon,NPP_Forest_by_flux2$lat,NPP_Forest_by_flux2$z), FUN=mean, na.rm=TRUE) #site-mean

#check
#analyse_modobs2(forest_site2,"pred_gpp", "GPP",type = "points")
ggplot(data=NPP_Forest_by_flux2, aes(x=new_gpp, y=GPP)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(GPP~new_gpp,NPP_Forest_by_flux2))

#analyse_modobs2(forest_site2,"pred_npp", "TNPP_1",type = "points")
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_npp, y=TNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(TNPP_1~pred_npp,NPP_Forest_by_flux2))

#analyse_modobs2(forest_site2,"pred_anpp", "ANPP_2",type = "points")
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_anpp, y=ANPP_2)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(ANPP_2~pred_anpp,NPP_Forest_by_flux2))

#analyse_modobs2(forest_site2,"pred_lnpp", "NPP.foliage",type = "points")
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_lnpp, y=NPP.foliage)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(NPP.foliage~pred_lnpp,NPP_Forest_by_flux2))

#analyse_modobs2(forest_site,"pred_wnpp", "NPP.wood",type = "points")
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_wnpp, y=NPP.wood)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(NPP.wood~pred_wnpp,NPP_Forest_by_flux2))

#analyse_modobs2(forest_site2,"pred_bnpp", "BNPP_1",type = "points")
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_bnpp, y=BNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(BNPP_1~pred_bnpp,NPP_Forest_by_flux2))

#analyse_modobs2(forest_site,"pred_lnf", "lnf_obs",type = "points") 
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_lnf, y=lnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(lnf_obs_final~pred_lnf,NPP_Forest_by_flux2))

#without second interpolation --> use this!
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_lnf, y=lnf_obs_org)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()+My_Theme
summary(lm(lnf_obs_org~pred_lnf,NPP_Forest_by_flux2))

#wnf - assuming constant wood/cn = 97
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_wnf, y=wnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(wnf_obs_final~pred_wnf,NPP_Forest_by_flux2))

#bnf - assuming constant root/cn = 122
ggplot(data=NPP_Forest_by_flux2, aes(x=pred_bnf, y=bnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(bnf_obs_final~pred_bnf,NPP_Forest_by_flux2))

#bnf - assuming constant root/cn = 122
NPP_Forest_by_flux2$pred_nuptake <- NPP_Forest_by_flux2$pred_lnf + NPP_Forest_by_flux2$pred_bnf + NPP_Forest_by_flux2$pred_wnf
NPP_Forest_by_flux2$obs_nuptake <- NPP_Forest_by_flux2$lnf_obs_org + NPP_Forest_by_flux2$bnf_obs_final + NPP_Forest_by_flux2$wnf_obs_final

ggplot(data=NPP_Forest_by_flux2, aes(x=pred_nuptake, y=obs_nuptake)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Prediction")+ylab("Observation")+theme_classic()
summary(lm(obs_nuptake~pred_nuptake,NPP_Forest_by_flux2))

save.image(file = "/Users/yunpeng/data/NPP_final/Forest_site_simulation.Rdata")

#now, adding Nuptake from additional source
NPP_Forest_by_flux2_Nuptake <- subset(NPP_Forest_by_flux2,obs_nuptake>0)
NPP_Forest_by_flux2_Nuptake <- NPP_Forest_by_flux2_Nuptake[,c("lon","lat","pred_nuptake","obs_nuptake")]
NPP_Forest_by_flux2_Nuptake$method <- "tissues C allocations and C/N"
Nuptake_new <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nuptake_gcme/All_Nuptake.csv")
Nuptake_new <- Nuptake_new[,c("lon","lat","pred_nuptake","obs_nuptake","method")]

Nuptake_final <- dplyr::bind_rows(NPP_Forest_by_flux2_Nuptake, Nuptake_new)

ggplot(data=Nuptake_final, aes(x=pred_nuptake, y=obs_nuptake)) +
  geom_point(aes(x=pred_nuptake, y=obs_nuptake,color=factor(method)))+geom_abline(intercept=0,slope=1)+geom_smooth(aes(x=pred_nuptake, y=obs_nuptake,color=factor(method)),method = "lm", se = TRUE)+
  xlab("Predicted N uptake")+ylab("Measured N uptake")+theme_classic()+My_Theme

ggplot(data=Nuptake_final, aes(x=pred_nuptake, y=obs_nuptake)) +
  geom_point(aes(x=pred_nuptake, y=obs_nuptake,color=factor(method)))+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Predicted N uptake")+ylab("Measured N uptake")+theme_classic()+My_Theme
summary(lm(obs_nuptake~pred_nuptake,Nuptake_final))

Nuptake_final_miner <- subset(Nuptake_final,method=="Net minerlization (Finzi paper)")
ggplot(data=Nuptake_final_miner, aes(x=pred_nuptake, y=obs_nuptake)) +
  geom_point(aes(x=pred_nuptake, y=obs_nuptake,color=factor(method)))+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Predicted N uptake")+ylab("Measured N uptake")+theme_classic()+My_Theme
summary(lm(obs_nuptake~pred_nuptake,Nuptake_final_miner))

ggplot(data=Nuptake_final, aes(x=pred_nuptake, y=obs_nuptake)) +
  geom_point(aes(x=pred_nuptake, y=obs_nuptake,color=factor(method)))+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  xlab("Predicted N uptake")+ylab("Measured N uptake")+theme_classic()+My_Theme
summary(lm(obs_nuptake~pred_nuptake,Nuptake_final))

newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)

#net minerlization
points(Nuptake_final$lon,Nuptake_final$lat, col="red", pch=16,cex=1)
points(Nuptake_final$lon[1:18],Nuptake_final$lat[1:18], col="green", pch=16,cex=1)
points(Nuptake_final$lon[259:346],Nuptake_final$lat[259:346], col="blue", pch=16,cex=1)
