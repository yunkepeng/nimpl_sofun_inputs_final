library(raster)
library(ncdf4)
library(dplyr)
library(maps)

#prepare a function for inputting monthly CRU data (tmn, tmx, vap and pre)
monthly_CRU <- function(variable,start_yr,end_yr,location){
  ncin <- nc_open(location)
  lon <- ncvar_get(ncin,"lon")
  nlon <- dim(lon)
  lat<-ncvar_get(ncin,"lat")
  nlat <- dim(lat)
  var_all <-ncvar_get(ncin,variable)
  nc_close(ncin)
  var_all2 <- as.vector(var_all)
  var_all3 <- matrix(var_all2, nrow = nlon * nlat, ncol = (2016-1901+1)*12) #as shown in file, 1392 months overall from 1901 to 2016 
  
  first_one <- (start_yr-1900)*12-11
  end_one <- (end_yr-1900)*12
  
  var_years <- var_all3[,c(first_one:end_one)]
  lonlat <- expand.grid(lon, lat)
  output_CRU <- cbind(lonlat,var_years)
  return(output_CRU)
}

monthly_tmn <- monthly_CRU("tmn",1980,2016,"/Volumes/My Passport/data/cru/ts_4.01/cru_ts4.01.1901.2016.tmn.dat.nc")
monthly_tmx <- monthly_CRU("tmx",1980,2016,"/Volumes/My Passport/data/cru/ts_4.01/cru_ts4.01.1901.2016.tmx.dat.nc")
monthly_vap <- monthly_CRU("vap",1980,2016,"/Volumes/My Passport/data/cru/ts_4.01/cru_ts4.01.1901.2016.vap.dat.nc")
monthly_pre <- monthly_CRU("pre",1980,2016,"/Volumes/My Passport/data/cru/ts_4.01/cru_ts4.01.1901.2016.pre.dat.nc")

#prepare a function for inputting monthly radiation data - based on mean of daily SWdown

#have a check about file - to see no trash files!
all_list <- list.files("/Volumes/My Passport/data/watch_wfdei/SWdown_daily/",full.names = T)
all_list[1] #197901
all_list[480] #201812

monthly_WFDEI <- function(variable,start_yr,end_yr,location){
  all_list <- list.files(location,full.names = T)
  first_one <- (start_yr-1978)*12-11
  end_one <- (end_yr-1978)*12
  years_list <- all_list[first_one:end_one]
  output_allyears <- data.frame(matrix(NA))
  #calculate average of each monthly data, within the selected years
  
  for (x in 1:length(years_list)){
    ncin <- nc_open(years_list[x])
    tstep<-ncvar_get(ncin,"timestp")
    n_days <- length(tstep)
    var_all <-ncvar_get(ncin,variable)
    nc_close(ncin)
    var_all2 <- as.vector(var_all)
    var_all3 <- matrix(var_all2, nrow = 720 * 360, ncol = n_days)
    var_all4 <- rowMeans(var_all3,na.rm = TRUE)
    output_allyears[1:259200,x] <- var_all4
  }
  #get coordinates info, and combine with final monthly average within selected years
  ncin <- nc_open(years_list[1])
  lon <- ncvar_get(ncin,"lon")
  lat<-ncvar_get(ncin,"lat")
  nc_close(ncin)
  
  lonlat <- expand.grid(lon, lat)
  output_WFDEI <- cbind(lonlat,output_allyears)
}

monthly_radi <- monthly_WFDEI("SWdown",1980,2016,"/Volumes/My Passport/data/watch_wfdei/SWdown_daily/")

save.image(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")



#below is what we copied to prediction.Rmd, to output climate prediction fields 
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")
library(rbeni)
library(dplyr)
library(raster)
library(ncdf4)
library(dplyr)
library(maps)
library(rgdal)

#calculate growth temperature
total_month <- (2016-1980+1) * 12

#solar declination from Jan to Dec: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/212421E8457A41C47B64D809BDEF53CA/9780511845727apx7_p351-352_CBO.pdf/solar_geometry_and_radiation_approximations.pdf
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

lat <- monthly_tmn[,2]
xx <- data.frame(matrix(, nrow=259200, ncol=total_month))
output_Tg <- data.frame(matrix(, nrow=259200, ncol=total_month))

for (a in 1:12){ 
  month_no <- seq(from = 1, to = total_month, by = 12)+a-1
  xx[1:259200,month_no]<- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)
#the percentage of tmx was dominated overall

output_Tg<-monthly_tmx[,3:(total_month+2)]*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ monthly_tmn[,3:(total_month+2)]*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
output_Tg[output_Tg =="NaN"] <- NA
output_Tg[output_Tg < 0] <- NA

output_vpd <- 0.611*exp(17.27*(output_Tg)/((output_Tg)+237.3))-monthly_vap[,3:(total_month+2)]*0.1 #vap in hPa
output_vpd[output_vpd =="NaN"] <- NA
output_vpd[output_vpd < 0] <- NA

output_PPFD <- monthly_radi[,3:(total_month+2)]*0.5*4.6 +output_Tg - output_Tg # here + Tg - Tg means it considered 
output_PPFD[output_PPFD =="NaN"] <- NA

#now, all done. Primarily calculated each year's average - then avaerage of years. Calculate na.rm=TRUE average.
Tg <- data.frame(matrix(, nrow=259200, ncol=12))
vpd <- data.frame(matrix(, nrow=259200, ncol=12))
PPFD <- data.frame(matrix(, nrow=259200, ncol=12))

#firstly calculate based on years (37 years to one year) - then calculate monthly average of one year
for (a in 1:12){ 
  month_no <- seq(from = 1, to = total_month, by = 12)+a-1
  Tg[1:259200,a]<- rowMeans(output_Tg[,month_no],na.rm = TRUE)
  vpd[1:259200,a]<- rowMeans(output_vpd[,month_no],na.rm = TRUE)
  PPFD[1:259200,a]<- rowMeans(output_PPFD[,month_no],na.rm = TRUE)
}

Tg <- rowMeans(Tg,na.rm = TRUE)
vpd <- rowMeans(vpd,na.rm = TRUE)
PPFD <- rowMeans(PPFD,na.rm = TRUE)

lonlat <- monthly_tmn[,1:2]

Tg_final <- cbind(lonlat,Tg)
names(Tg_final) <- c("lon","lat","Tg")
Tg_final$Tg[Tg_final$Tg=="NaN"] <- NA
Tg_df_output <- Tg_final

vpd_final <- cbind(lonlat,vpd)
names(vpd_final) <- c("lon","lat","vpd")
vpd_final$vpd[vpd_final$vpd=="NaN"] <- NA
vpd_df_output <- vpd_final

PPFD_final <- cbind(lonlat,PPFD)
names(PPFD_final) <- c("lon","lat","PPFD")
PPFD_final$PPFD[PPFD_final$PPFD=="NaN"] <- NA
PPFD_df_output <- PPFD_final

#prepare lon and lat
library(ncdf4)
ncin <- nc_open("~/data/watch_wfdei/WFDEI-elevation.nc")
lon <- ncvar_get(ncin,"lon")
lat<-ncvar_get(ncin,"lat")

#output PPFD nc file - In Euler its path is same: "~/data/nimpl_sofun_inputs/map/Final_ncfile"
PPFD_nc <- list(df_to_grid(PPFD_df_output,varnam = "PPFD", lonnam = "lon", latnam = "lat"))
names(PPFD_nc) <- "PPFD"
varams = "PPFD"
test <- list(lon,lat,PPFD_nc,varams)
names(test) <- c("lon","lat","vars","varams")
write_nc2(test,varnams = "PPFD",long_name = "PPFD",units = "umol/m2/s",
          path = "~/data/nimpl_sofun_inputs/map/Final_ncfile/PPFD.nc")

#output Tg nc file - In Euler its path is same: "~/data/nimpl_sofun_inputs/map/Final_ncfile"
Tg_nc <- list(df_to_grid(Tg_df_output,varnam = "Tg", lonnam = "lon", latnam = "lat"))
names(Tg_nc) <- "Tg"
varams = "Tg"
test <- list(lon,lat,Tg_nc,varams)
names(test) <- c("lon","lat","vars","varams")
write_nc2(test,varnams = "Tg",long_name = "Growth temperature",units = "Degree Celcius",
          path = "~/data/nimpl_sofun_inputs/map/Final_ncfile/Tg.nc")

#output vpd nc file - In Euler its path is same: "~/data/nimpl_sofun_inputs/map/Final_ncfile"
vpd_nc <- list(df_to_grid(vpd_df_output,varnam = "vpd", lonnam = "lon", latnam = "lat"))
names(vpd_nc) <- "vpd"
varams = "vpd"
test <- list(lon,lat,vpd_nc,varams)
names(test) <- c("lon","lat","vars","varams")
write_nc2(test,varnams = "vpd",long_name = "vapor prssure deficient",units = "KPa",
          path = "~/data/nimpl_sofun_inputs/map/Final_ncfile/vpd.nc")



library(raster)

coordinates(PPFD_final) <- ~lon+lat 
gridded(PPFD_final) <- TRUE
r3 <- raster(PPFD_final, "PPFD") 
plot(r3)

coordinates(Tg_final) <- ~lon+lat 
gridded(Tg_final) <- TRUE
r3 <- raster(Tg_final, "Tg") 
plot(r3)

coordinates(vpd_final) <- ~lon+lat 
gridded(vpd_final) <- TRUE
r3 <- raster(vpd_final, "vpd") 
plot(r3)
