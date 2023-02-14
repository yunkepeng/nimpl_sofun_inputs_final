#prepare sites
rm(list=ls())
devtools::load_all("/Users/yunpeng/yunkepeng/latest_packages/rbeni/") 
library(rworldmap)
library(spgwr)
library(readr)
#input site info
allsites <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Nmin_dataset.csv")
gwr_sites <- aggregate(allsites,by=list(allsites$lon,allsites$lat,allsites$z,allsites$Begin_year,allsites$End_year), FUN=mean, na.rm=TRUE)
gwr_sites <- gwr_sites[,c("lon","lat","z","Begin_year","End_year")]

gwr_sites$year_start <- gwr_sites$Begin_year
gwr_sites$year_end <- gwr_sites$End_year
gwr_sites$year_start[gwr_sites$Begin_year<=1980] <- 1980
gwr_sites$year_end[gwr_sites$End_year<=1980] <- 1989

dim(gwr_sites)
summary(gwr_sites)

#input tmx, tmn, vap, pre and radi
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")

total_month <- (2016-1980+1) * 12

elev <- as.data.frame(nc_to_df(read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc"), varnam = "elevation"))

monthly_tmn$lon <- elev$lon
monthly_tmn$lat <- elev$lat
monthly_tmn$z <- elev$elevation
monthly_tmn <- monthly_tmn[,c(3: (total_month+5))] #lon, lat, z in last 3
monthly_tmn <- subset(monthly_tmn,z>=0)

#just check if points are all correct, for available lon+lat+z after cbind climates data - we disregard those grids with z=NA or negative values (just sea..), because we input elevation necessarily in gwr below! 
library(maps)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(monthly_tmn$lon,monthly_tmn$lat, col="red", pch=16,cex=1)
#ok

monthly_tmx$lon <-  elev$lon
monthly_tmx$lat <- elev$lat
monthly_tmx$z <- elev$elevation
monthly_tmx <- monthly_tmx[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_tmx <- subset(monthly_tmx,z>=0)

monthly_vap$lon <-  elev$lon
monthly_vap$lat <- elev$lat
monthly_vap$z <- elev$elevation
monthly_vap <- monthly_vap[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_vap <- subset(monthly_vap,z>=0)

monthly_pre$lon <-  elev$lon
monthly_pre$lat <-  elev$lat
monthly_pre$z <- elev$elevation
monthly_pre <- monthly_pre[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_pre <- subset(monthly_pre,z>=0)

monthly_radi$lon <- elev$lon
monthly_radi$lat <- elev$lat
monthly_radi$z <- elev$elevation
monthly_radi <- monthly_radi[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_radi <- subset(monthly_radi,z>=0)

#alpha extracted from directly from monthly map (as used in prediction field)
alphalist <- list.files(path = "~/data/alpha/data_orig/",full.names = T) # here this file includes 116 data - that is the annual alpha from 1901 to 2016
empty_alpha <- data.frame(matrix(NA)) 
#calculate alpha within 37 years (1980-2016)
alphalist[80]
alphalist[116]
for (i in 80:116){ #here 82:111 means 1980-2016
  load(file = alphalist[i])
  empty_alpha[1:259200, (((i-79)*12-11):((i-79)*12))] <- SP_result_monthly
}

empty_alpha$lon <- elev$lon
empty_alpha$lat <- elev$lat
empty_alpha$z <- elev$elevation

empty_alpha <- subset(empty_alpha,z>=0)
dim(empty_alpha) #lon, lat, z in last 3

d <- 1.5 # range of surrounded grids (in degree) when applying each gwr

#start gwr 
library(spgwr)
library(raster)

#create a function here
gwr_methods <- function(sites,monthly_grid){
  for (i in c(1:nrow(sites))){
    print(i)
    start_yr <- sites$year_start[i]
    end_yr <- sites$year_end[i]
    lon_site <- sites$lon[i]
    lat_site <- sites$lat[i]
    z_site <- sites$z[i]
    
    climate_distance <- subset(monthly_grid,lon>(lon_site-d)&lon<(lon_site+d)&lat>(lat_site-d)&lat<(lat_site+d))
    climate_distance_yr <- climate_distance[,((start_yr-1979)*12-11):((end_yr-1979)*12)]
    number_of_month <- ncol(climate_distance_yr)
    
    climate_distance_yr_xyz <- cbind(climate_distance[,c("lon","lat","z")],climate_distance_yr)
    # re-create column name for multiple month, which is important when applying gwr within a for loop
    names(climate_distance_yr_xyz) <- c("lon","lat","z",paste("a", 1:number_of_month, sep=""))
    
    #filter the condition when without available climates data after subset (filter within a distance of 1.5)
    if (nrow(climate_distance) == 0){
      sites[i,8:(7+number_of_month)] <- NA
    } else {
      sites[i,8:(7+number_of_month)] <- NA
      coordinates(climate_distance_yr_xyz) <- c("lon","lat")
      gridded(climate_distance_yr_xyz) <- TRUE
      
      input_sites <- sites[i,c("lon","lat","z")]
      coordinates(input_sites) <- c("lon","lat")
      for (a in c(1:number_of_month)){
        #a +5 means, add outputted monthly tmn (within measurement year) directly after lon+lat+z+Begin_year+End_year.
        sites[i,a+7] <- (gwr( get(paste("a",a, sep="")) ~ z, climate_distance_yr_xyz, bandwidth = 1.06, fit.points = input_sites,predictions=TRUE))$SDF$pred
      }
    }
  }
  output <- sites
  return(output)
}


alpha_output <- gwr_methods(gwr_sites,empty_alpha)

tmx_output <- gwr_methods(gwr_sites,monthly_tmx)

tmn_output <- gwr_methods(gwr_sites,monthly_tmn)

vap_output <- gwr_methods(gwr_sites,monthly_vap)

pre_output <- gwr_methods(gwr_sites,monthly_pre)

radi_output <- gwr_methods(gwr_sites,monthly_radi)

#output it temporaily
csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_site_newphy.csv")
write_csv(gwr_sites, path = csvfile)

csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_alpha_newphy.csv")
write_csv(alpha_output, path = csvfile)

csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_tmx_newphy.csv")
write_csv(tmx_output, path = csvfile)

csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_tmn_newphy.csv")
write_csv(tmn_output, path = csvfile)

csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_vap_newphy.csv")
write_csv(vap_output, path = csvfile)

csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_pre_newphy.csv")
write_csv(pre_output, path = csvfile)

csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_radi_newphy.csv")
write_csv(radi_output, path = csvfile)

# now, calculate Tg, vpd and PPFD 
dim(radi_output)
tmx_site <- tmx_output[,8:247]
tmn_site <- tmn_output[,8:247]
vap_site <- vap_output[,8:247]
pre_site <- pre_output[,8:247]
radi_site <- radi_output[,8:247]
alpha_site <- alpha_output[,8:247]

#1. Tg
#solar declination from Jan to Dec
s1 <- (-23.1+ -17.3)/2
s2 <- (-17.3 + -8)/2
s3 <- (-8 + 4.1)/2
s4 <- (4.1 + 14.8)/2
s5 <- (14.8 + 21.9)/2
s6 <- (21.9 + 23.2)/2
s7 <- (23.2 + 18.3)/2
s8 <- (18.3 + 8.6)/2
s9 <-  (8.6 + -2.8)/2
s10 <- (-2.8 + -14.1)/2
s11 <- (-14.1 + -21.6)/2
s12 <- (-21.6 + -23.1)/2

s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
lat <- gwr_sites$lat
xx <- data.frame(matrix(nrow=nrow(gwr_sites), ncol=240))
output_Tg <- data.frame(matrix(nrow=nrow(gwr_sites), ncol=240))
#xx = acos(h), h = hour angle of the sun
for (a in 1:12){ 
  month_no <- seq(from = 1, to = 240, by = 12)+a-1
  xx[1:nrow(gwr_sites),month_no]<- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)

#the percentage of tmx was dominated overall
Tg_site <- tmx_site*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ tmn_site*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
Tg_site[Tg_site =="NaN"] <- NA
Tg_site[Tg_site < 0] <- NA

vpd_site <- 0.611*exp(17.27*(Tg_site)/((Tg_site)+237.3))-vap_site*0.1 #vap in hPa
vpd_site[vpd_site =="NaN"] <- NA
vpd_site[vpd_site < 0] <- NA

PPFD_site <- radi_site*0.5*4.6 + Tg_site - Tg_site # here + Tg - Tg means it considered growth season
PPFD_site[PPFD_site =="NaN"] <- NA
PPFD_site[PPFD_site < 0] <- NA

alpha_Tg_site <- alpha_site+ Tg_site - Tg_site # here + Tg - Tg means it considered growth season
alpha_Tg_site[alpha_Tg_site >1] <- 1
alpha_Tg_site[alpha_Tg_site < 0] <- NA


#merged now
gwr_sites$alpha_sites <- rowMeans(alpha_Tg_site,na.rm=TRUE)
gwr_sites$PPFD_sites <- rowMeans(PPFD_site,na.rm=TRUE)
gwr_sites$Tg_sites <- rowMeans(Tg_site,na.rm=TRUE)
gwr_sites$vpd_sites <- rowMeans(vpd_site,na.rm=TRUE)

summary(gwr_sites)
csvfile <- paste("/Users/yunpeng/data/NPP_Yunke/predictors/climates_sites.csv")
write_csv(gwr_sites[,c("lon","lat","z","Begin_year","End_year","year_start","year_end","alpha_sites","PPFD_sites","Tg_sites","vpd_sites")], path = csvfile)
