#create a function to pre-processing age, so that output each bins of age (from 4-pft) to file in /Users/yunpeng/data/nimpl_sofun_inputs/map/Final_ncfile
library(raster)
library(ncdf4)
library(dplyr)
library(maps)
library(rgdal)
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")

#set the path and filename
ncfname <- paste ("~/data/GFAD/data_orig/", "GFAD_V1-1", ".nc", sep="")
#ncfname <- paste ("E:/C-N cycling/Carbon allocation/672 sites analysis/stand age/GFAD2/", "GFAD_V1-1", ".nc", sep="")

#1. open a netCDF file and check variable
ncin <- nc_open(ncfname)

#2. Get coordinate and value
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon) 
lat<-ncvar_get(ncin,"lat")
nlat<-dim(lat)
age <-ncvar_get(ncin,"age")
nc_close(ncin)

pre.vec.long <- as.vector(age)
pre.mat <- matrix(pre.vec.long, nrow = nlon * nlat, ncol = 4*15) # 259200 columns * 60 rows - here 60 represents 4 PFTs and 15 age classes.

lonlat <- expand.grid(lon, lat)
lonlat$age_5 <- rowSums(pre.mat[,1:4])
lonlat$age_15 <- rowSums(pre.mat[,5:8])
lonlat$age_25 <- rowSums(pre.mat[,9:12])
lonlat$age_35 <- rowSums(pre.mat[,13:16])
lonlat$age_45 <- rowSums(pre.mat[,17:20])
lonlat$age_55 <- rowSums(pre.mat[,21:24])
lonlat$age_65 <- rowSums(pre.mat[,25:28])
lonlat$age_75 <- rowSums(pre.mat[,29:32])
lonlat$age_85 <- rowSums(pre.mat[,33:36])
lonlat$age_95 <- rowSums(pre.mat[,37:40])
lonlat$age_105 <- rowSums(pre.mat[,41:44])
lonlat$age_115 <- rowSums(pre.mat[,45:48])
lonlat$age_125 <- rowSums(pre.mat[,49:52])
lonlat$age_135 <- rowSums(pre.mat[,53:56])
lonlat$age_145 <- rowSums(pre.mat[,57:60])

age_input <- as.data.frame(lonlat)
age_input$lon <- age_input$Var1
age_input$lat <- age_input$Var2
age_output <- function(df,names){
  age_input <- age_input[,c("lon","lat",names)]
  names(age_input) <- c("lon","lat","age")
  coordinates(age_input) <- ~lon+lat 
  gridded(age_input) <- TRUE
  r <- raster(age_input,"age")
  
  ncfname <- paste ("~/data/landmasks/", "srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant", ".nc", sep="")
  #ncfname <- paste ("D:/PhD/nimpl_sofun_inputs/Data/srex/", "srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant", ".nc", sep="")
  ncin <- nc_open(ncfname)
  
  lon <- ncvar_get(ncin,"lon")
  nlon <- dim(lon) 
  lat<-ncvar_get(ncin,"lat")
  nlat <- dim(lat)
  area <- ncvar_get(ncin,"srex_mask")
  nc_close(ncin)
  
  area1 <- as.vector(area)
  
  lonlat <- expand.grid(lon, lat)
  continent_input <- as.data.frame(cbind(lonlat,area1))
  names(continent_input) <- c("lon","lat","area")
  
  #merge continents and age
  continent_age <- merge(continent_input, age_input, by=c("lon","lat"), all.x = T, sort=F )
  summary(continent_age)
  
  #3. Input fAPAR, and merge it with continent_age, and convert fAPAR's NA grid to those age = NA
  ncfname <- paste ("~/data/fAPAR/fAPAR3g_v2/", "fAPAR3g_v2_1982_2016_FILLED", ".nc", sep="")
  #ncfname <- paste ("D:/PhD/nimpl_sofun_inputs/Data/fAPAR/", "fAPAR3g_v2_1982_2016_FILLED", ".nc", sep="")
  dname <- "FAPAR_FILLED"
  
  ncin <- nc_open(ncfname)
  
  lon <- ncvar_get(ncin,"LON")
  nlon <- dim(lon) 
  
  lat<-ncvar_get(ncin,"LAT")
  nlat <- dim(lat)
  
  FAPAR <- ncvar_get(ncin,"FAPAR_FILLED")
  nc_close(ncin)
  
  pre.vec.long <- as.vector(FAPAR)
  
  pre.mat <- matrix(pre.vec.long, nrow = nlon * nlat, ncol = 420) # 259200 * 420 (420 is 35 years: 1982-2016)
  
  fAPAR1982_2011 <- pre.mat[,c(1:360)] #extrat fAPAR from 1982-2011
  
  final_fAPAR <- rowMeans(fAPAR1982_2011,na.rm = TRUE) # calculate means
  
  lonlat <- expand.grid(lon, lat)
  fAPAR_input <- as.data.frame(cbind(lonlat,final_fAPAR))
  names(fAPAR_input) <- c("lon","lat","fAPAR")

  #here we merge fAPAR with age and continent data.
  continent_age_fAPAR <- merge(continent_age, fAPAR_input, by=c("lon","lat"), all.x = T, sort=F )

  #when fAPAR = NA, then those grid's age also = NA. After this, age = 0 only had one possibility, that is, missing data. We will subset them later on and replace them to mean value of local continent.
  for (i in 1:nrow(continent_age_fAPAR)){
    if (is.na(continent_age_fAPAR$fAPAR[i]) == TRUE){
      continent_age_fAPAR$new_age[i]<- NA} else { 
        continent_age_fAPAR$new_age[i]<- continent_age_fAPAR$age[i]} 
  }
  
  continent_age_fAPAR2 <- continent_age_fAPAR[,c("lon","lat","new_age")]
  
  coordinates(continent_age_fAPAR2) <- ~lon+lat 
  gridded(continent_age_fAPAR2) <- TRUE
  r <- raster(continent_age_fAPAR2, "new_age") 

  #4. Now, start interpolating each continent's missing data (when age = 0) as mean value of local continent
  age_final <- as.data.frame(continent_age_fAPAR[,c("lon","lat","area","new_age")])
  
  mylist <- vector(mode = "list", length =25) # for specifying gridded data at certain degrees
  
  for (i in c(1:8,10:24)){ # area 9 and 25 did not have such problem (missing data of stand-age; i.e. age = 0)
    empty_age <- subset(age_final, (area == i & new_age == 0)) #find missing grids (new_age ==0, classfied as missing data)
    available_age <- subset(age_final, (area == i & new_age > 0)) #find non-missing grids and calculate mean
    empty_age$final_age <- mean(available_age$new_age)
    mylist[[i]] <- empty_age
    #print(i)
  }
  
  all_age <- as.data.frame(do.call(rbind, mylist))

  #South AUS (area = 26), however, was interpolated by area =25 's mean stand-age value
  
  South_AUS_empty_age <- subset(age_final, (area == 26 & new_age == 0) )
  North_AUS_available_age <- subset(age_final, (area == 25 & new_age > 0))
  South_AUS_empty_age$final_age <- mean(North_AUS_available_age$new_age)
  
  all_age2 <- rbind(all_age,South_AUS_empty_age)
  
  all_age3 <- all_age2[,c("lon","lat","final_age")]
  
  #5. merge new interpolated subset of age dataframe into last complete dataframe, and generating findal data to be used.
  final <- merge(continent_age_fAPAR, all_age3, by=c("lon","lat"), all.x = T,sort = T)
  
  for (i in 1:nrow(final)){
    if (is.na(final$final_age[i]) == TRUE){
      final$age_used[i] <- final$new_age[i]} else {
        final$age_used[i] <- final$final_age[i]} 
  }
  
  final2 <- final[,c("lon","lat","age_used")]
  names(final2) <- c("lon","lat","age")
  
  final_data <- final2
  
  coordinates(final2) <- ~lon+lat 
  gridded(final2) <- TRUE
  r <- raster(final2, "age") 

  #6. convert dataframe to nc file as input in nimpl project
  outlier <- (subset(final_data,age ==0))
  
  final_data$age[final_data$age == 0] <- NA
  #now final_data is ready to be used, let's firstly change its lat, lon order correctly, then convert all NA to 9999
  
  final_data2 <- final_data[order(final_data[,2],final_data[,1]),]
  
  df_age <- final_data2
  
  #output age map
  coordinates(final_data2) <- ~lon+lat 
  gridded(final_data2) <- TRUE
  rage <- raster(final_data2, "age") 
  plot(rage)
  
  #output nc file - In Euler its path is same: "~/data/nimpl_sofun_inputs/map/Final_ncfile"
  summary(df_age)
  age_nc <- list(df_to_grid(df_age,varnam = "age", lonnam = "lon", latnam = "lat"))
  names(age_nc) <- "age"
  varams = "age"
  test <- list(lon,lat,age_nc,varams)
  names(test) <- c("lon","lat","vars","varams")
  write_nc2(test,varnams = "age",long_name = "stand age",units = "years",
            path = paste("~/data/nimpl_sofun_inputs/map/Final_ncfile/",names,".nc",sep=""))
  return(df_age)
}
names(age_input)

age_5 <- age_output(age_input,"age_5")
age_15 <- age_output(age_input,"age_15")
age_25 <- age_output(age_input,"age_25")
age_35 <- age_output(age_input,"age_35")
age_45 <- age_output(age_input,"age_45")
age_55 <- age_output(age_input,"age_55")
age_65 <- age_output(age_input,"age_65")
age_75 <- age_output(age_input,"age_75")
age_85 <- age_output(age_input,"age_85")
age_95 <- age_output(age_input,"age_95")
age_105 <- age_output(age_input,"age_105")
age_115 <- age_output(age_input,"age_115")
age_125 <- age_output(age_input,"age_125")
age_135 <- age_output(age_input,"age_135")
age_145 <- age_output(age_input,"age_145")

# a test to show wrong age
a1 <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/age_105.nc"),
  varnam = "age"))

plot_map3(na.omit(a1[,c("lon","lat","myvar")]),
          varnam = "myvar",latmin = -65, latmax = 85)
