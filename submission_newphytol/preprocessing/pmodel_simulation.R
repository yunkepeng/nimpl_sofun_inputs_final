rm(list=ls())

#load(file = "/Users/yunpeng/data/NPP_final/Forest_site_simulation.Rdata")
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")

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
#NPP_Forest <- read.csv("/Users/yunpeng/data/NPP_final/fpar_name/forest_fpar_name.csv")
NPP_Forest <- read.csv("/Users/yunpeng/data/NPP_Yunke/simulated_gpp/forest_fpar_name.csv")

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

#we have newly interpolate their fapar primarily based on n_focal = 1, then n_focal = 2, and saved it in "/Users/yunpeng/data/forest_npp/reprocessing_fpar_raw/"
#the code of this is available at L90-110 in forest/Reprocessing_fpar_climates_forest.R

#now, reprocessing such values - by updating such fapar ###this code (from 79 to 86) add addtional data from n_focal = 1 or 2 (if not using them the don't run!)
#fapar_df_new <- list.files("~/data/forest_npp/reprocessing_fpar_raw/",full.names = T)
fapar_df_new <- list.files("~/data/NPP_final/reprocessing_fpar_raw2/",full.names = T)

for (i in 1:(length(fapar_df_new)-1)){
  df1 <- read.csv(fapar_df_new[i])
  df1$date <- as.Date(df1$date)
  df1 <- df1[!(format(df1$date,"%m") == "02" & format(df1$date, "%d") == "29"), , drop = FALSE]
  df2 <- df1[,c("date","modisvar_filled")]
  assign(substr(sub('.*daily_', '', fapar_df_new[i]),1,nchar(sub('.*daily_', '', fapar_df_new[i]))-4), df2) 
}

#check climate forcing missing data
empty_vec <- c()

#check existed climate files
for (i in 1:(length(forcing_df))){
  empty_vec[i] <- as.numeric(gsub("[^0-9]", "",  forcing_df[i]))
}

diff <- setdiff(1:935, empty_vec)
diff
#NPP_F338,NPP_390

#totally 5 sites were missing:
all_na_points <- c(na_fapar$sitename,c("NPP_F338","NPP_390"))
#NPP_F556, NPP_F697, NPP_F700 (due to fapar orig missing in ingestr_fpar261,263,264) and NPP_F338, NPP_F390 (due to climate forcing missing)

#2. forcing - combing fapar and climates into a df.
#there are some extra csv from forcing files (e.g. Sara Vicca_flux), which was not used in the final version of rsofun -because such additional data in 47 points don't have flux data, but only have biomass data (not useful here).
#therefore, it would strictly follows data that it included in NPP_Forest
for (i in 1:length(forcing_df)){
  tryCatch({
    df1 <- read.csv(forcing_df[i])
    df1$date <- as.Date(df1$date)
    #below is important - if consistent then continue - if not then it may be from Sara Vicca_flux therefore not available
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
    start_date <- as.numeric(substr(forcing$date[1],1,4)) # check if it is consistent with org data
    end_date <- as.numeric(substr(forcing$date[length(forcing$date)],1,4))  # check if it is consistent with org data
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
    NPP_Forest[i,c("start_date")] <- start_date
    NPP_Forest[i,c("end_date")] <- end_date
  }, error=function(e){})} 

#check if time is consistent
summary(NPP_Forest$year_start -NPP_Forest$start_date)
summary(NPP_Forest$year_end -NPP_Forest$end_date)
subset(NPP_Forest,year_end!=end_date) #ok, but within range

#check missing gpp
subset(NPP_Forest,is.na(pred_gpp_c3)==TRUE)$sitename_fpar
length(subset(NPP_Forest,is.na(pred_gpp_c3)==TRUE)$sitename_fpar)
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(subset(NPP_Forest,is.na(pred_gpp_c3)==TRUE)$lon,subset(NPP_Forest,is.na(pred_gpp_c3)==TRUE)$lat, col="red", pch=16,cex=1)
#these points were missing, either due to fapar or climate forcing missing

#combine with current dataset
measurement <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Nmin_dataset.csv")
measurement_forest <- subset(measurement,pft=="Forest" &is.na(Nmin)==TRUE)
measurement_forest$year_start <- measurement_forest$Begin_year
measurement_forest$year_end <- measurement_forest$End_year

measurement_forest$year_start[measurement_forest$Begin_year<=1980] <- 1980
measurement_forest$year_end[measurement_forest$End_year<=1980] <- 1989
dim(measurement_forest)

merged_data <- NPP_Forest[,c("lon","lat","elv","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")]
merged_data_forest <- aggregate(merged_data,by=list(merged_data$lon,merged_data$lat,merged_data$elv,merged_data$year_start,merged_data$year_end), FUN=mean, na.rm=TRUE)[,c("lon","lat","elv","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")]
names(merged_data_forest) <- c("lon","lat","z","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")

#merged with existing measurements - all good
measurement_forest2 <- merge(measurement_forest,merged_data_forest,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
summary(measurement_forest2)

ggplot(data=subset(measurement_forest2,rep=="not_repeated"), aes(x=pred_gpp_c3, y=GPP)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  labs(y = ~paste("Forest ", GPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", GPP[pred.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,5000)+ylim(0,5000)


#grassland
#siteinfo_final <- as.data.frame(read_csv("/Users/yunpeng/data/NPP_final/fpar_name/grassland_fpar_name.csv"))
siteinfo_final <- as.data.frame(read_csv("/Users/yunpeng/data/NPP_Yunke/simulated_gpp/grassland_fpar_name.csv"))

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
siteinfo_final$pred_gpp_c3 <- NA
siteinfo_final$pred_gpp_c4 <- NA
siteinfo_final$max_vcmax25_c3 <- NA
siteinfo_final$max_vcmax25_c4 <- NA

#now, gpp-c3, gpp-c4, vcmax25max-c3, vcmax25-c4
for (i in 1:nrow(siteinfo_final)) {
  tryCatch({
    #c3 gpp
    forcing <- (eval(parse(text=(paste("final",siteinfo_final$sitename[i],sep="_")))))
    start_date <- as.numeric(substr(forcing$date[1],1,4)) # check if it is consistent with org data
    end_date <- as.numeric(substr(forcing$date[length(forcing$date)],1,4))  # check if it is consistent with org data
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
    
    siteinfo_final[i,c("start_date")] <- start_date
    siteinfo_final[i,c("end_date")] <- end_date
    
  }, error=function(e){})} 

#check if time is consistent - good
summary(siteinfo_final$year_start -siteinfo_final$start_date)
summary(siteinfo_final$year_end -siteinfo_final$end_date)

#combine with current dataset
measurement <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Nmin_dataset.csv")
measurement_grassland <- subset(measurement,pft!="Forest" &is.na(Nmin)==TRUE)
measurement_grassland$year_start <- measurement_grassland$Begin_year
measurement_grassland$year_end <- measurement_grassland$End_year

measurement_grassland$year_start[measurement_grassland$Begin_year<=1980] <- 1980
measurement_grassland$year_end[measurement_grassland$End_year<=1980] <- 1989

merged_grassland <- siteinfo_final[,c("lon","lat","elv","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")]
merged_data_grassland <- aggregate(merged_grassland,by=list(merged_grassland$lon,merged_grassland$lat,merged_grassland$elv,merged_grassland$year_start,merged_grassland$year_end), FUN=mean, na.rm=TRUE)[,c("lon","lat","elv","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")]
names(merged_data_grassland) <- c("lon","lat","z","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")

#merged with existing measurements - all good
measurement_grassland2 <- merge(measurement_grassland,merged_data_grassland,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)

ggplot(data=subset(measurement_grassland2,rep=="not_repeated"), aes(x=pred_gpp_c3, y=GPP)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  labs(y = ~paste("Grassland ", GPP[obs.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland ", GPP[pred.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,5000)+ylim(0,5000)

fill_missing <- unique(subset(measurement_grassland2,is.na(max_vcmax25_c3)==TRUE)[,c("lon","lat","z","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")])

#fill some missing data (some lacked because one measurement year (in old simulation converted to) show 1980-1989 while measurement data show 1980-1982, so that failed to merge)
fill_missing$pred_gpp_c3[fill_missing$lon==11.8] <-merged_grassland$pred_gpp_c3[merged_grassland$lon==11.8];fill_missing$max_vcmax25_c3[fill_missing$lon==11.8] <-merged_grassland$max_vcmax25_c3[merged_grassland$lon==11.8];
fill_missing$pred_gpp_c3[fill_missing$lon==36.50000&fill_missing$year_start==1980] <-merged_grassland$pred_gpp_c3[merged_grassland$lon==36.50000&merged_grassland$year_start==1980];fill_missing$max_vcmax25_c3[fill_missing$lon==36.50000&fill_missing$year_start==1980] <-merged_grassland$max_vcmax25_c3[merged_grassland$lon==36.50000&merged_grassland$year_start==1980]
fill_missing$pred_gpp_c3[fill_missing$lon==60.08000&fill_missing$year_start==1980] <-merged_grassland$pred_gpp_c3[merged_grassland$lon==60.08000&merged_grassland$year_start==1980];fill_missing$max_vcmax25_c3[fill_missing$lon==60.08000&fill_missing$year_start==1980] <-merged_grassland$max_vcmax25_c3[merged_grassland$lon==60.08000&merged_grassland$year_start==1980]
fill_missing$pred_gpp_c3[fill_missing$lon==116.55000&fill_missing$year_start==1980] <-merged_grassland$pred_gpp_c3[merged_grassland$lon==116.55000&merged_grassland$year_start==1980];fill_missing$max_vcmax25_c3[fill_missing$lon==116.55000&fill_missing$year_start==1980] <-merged_grassland$max_vcmax25_c3[merged_grassland$lon==116.55000&merged_grassland$year_start==1980]
fill_missing$pred_gpp_c3[fill_missing$lon==116.66000&fill_missing$year_start==1980] <-merged_grassland$pred_gpp_c3[merged_grassland$lon==116.66000&merged_grassland$year_start==1980];fill_missing$max_vcmax25_c3[fill_missing$lon==116.66000&fill_missing$year_start==1980] <-merged_grassland$max_vcmax25_c3[merged_grassland$lon==116.66000&merged_grassland$year_start==1980]
fill_missing$pred_gpp_c3[fill_missing$lon==116.74000&fill_missing$year_start==1980] <-merged_grassland$pred_gpp_c3[merged_grassland$lon==116.74000&merged_grassland$year_start==1980];fill_missing$max_vcmax25_c3[fill_missing$lon==116.74000&fill_missing$year_start==1980] <-merged_grassland$max_vcmax25_c3[merged_grassland$lon==116.74000&merged_grassland$year_start==1980]

#some measured year were missed to be collected so it has not merged sucessfuly - just fill it - basing on lon, lat, z only
merged_data_grassland2 <- na.omit(fill_missing)

#finally, site simulation of N uptake
#input Filzi's climates
Finzi <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nmin_Finzi/Nmin_Finzi.csv")
names(Finzi)[names(Finzi) == "Lat"] <- "lat"
names(Finzi)[names(Finzi) == "Long"] <- "lon"
devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/ingestr/")

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

Nmin_final <- Finzi_all_forest

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

merged_data_Nmin <- na.omit(unique(Nmin_final[,c("lon","lat","elv","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")]))
names(merged_data_Nmin) <- c("lon","lat","z","year_start","year_end","pred_gpp_c3","max_vcmax25_c3")
all_forcing <- dplyr::bind_rows(merged_data_forest, merged_data_grassland,merged_data_grassland2,merged_data_Nmin) 
summary(all_forcing)

csvfile <- paste("/Users/yunpeng/data/NPP_Yunke/simulated_gpp/site_simulated_gpp_vcmax.csv")
write_csv(all_forcing, path = csvfile)