library(raster)
library(ncdf4)
library(dplyr)
library(maps)
library(rgdal)
devtools::load_all("/Users/yunpeng/yunkepeng/latest_packages/rbeni/") 

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
lonlat$age <- rowSums(pre.mat[,1:4])*5+rowSums(pre.mat[,5:8])*15+rowSums(pre.mat[,9:12])*25+rowSums(pre.mat[,13:16])*35+rowSums(pre.mat[,17:20])*45+rowSums(pre.mat[,21:24])*55+rowSums(pre.mat[,25:28])*65+rowSums(pre.mat[,29:32])*75+rowSums(pre.mat[,33:36])*85+rowSums(pre.mat[,37:40])*95+rowSums(pre.mat[,41:44])*105+rowSums(pre.mat[,45:48])*115+rowSums(pre.mat[,49:52])*125+rowSums(pre.mat[,53:56])*135+rowSums(pre.mat[,57:60])*145
names(lonlat) <- c("lon","lat","age")

age_input <- as.data.frame(lonlat)
summary(age_input)
plot_map3(age_input, 
          varnam = "age",plot_title = "age",
          latmin = -65, latmax = 85)
dim(age_input)

plot_map3
# now it shows an empty grid in lower AUS, also in some other regions. We need to interpolate them based on local mean value.

#coordinates(age_input) <- ~lon+lat 
#gridded(age_input) <- TRUE
#r <- raster(age_input, "age") 
#plot(r) 


#2. Input continent map
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
summary(fAPAR_input)

#here we merge fAPAR with age and continent data.
continent_age_fAPAR <- merge(continent_age, fAPAR_input, by=c("lon","lat"), all.x = T, sort=F )
summary(continent_age_fAPAR)
dim(continent_age_fAPAR)

#when fAPAR = NA, then those grid's age also = NA. After this, age = 0 only had one possibility, that is, missing data. We will subset them later on and replace them to mean value of local continent.
continent_age_fAPAR$new_age <- NA
continent_age_fAPAR$new_age[is.na(continent_age_fAPAR$fAPAR)==FALSE] <- continent_age_fAPAR$age[is.na(continent_age_fAPAR$fAPAR)==FALSE]
continent_age_fAPAR$new_age[is.na(continent_age_fAPAR$fAPAR)==TRUE] <- NA

continent_age_fAPAR2 <- continent_age_fAPAR[,c("lon","lat","new_age")]


plot_map3(continent_age_fAPAR2, 
          varnam = "new_age",plot_title = "age before filled by AUS",
          latmin = -65, latmax = 85)

#we will start working from this map! Because it actually clarfies which grid stand-age = NA (in white), which grid stand-age has missed data (which temporaily considered = 0, and shown in grey, will be interpolated soon) .


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
dim(all_age)

#South AUS (area = 26), however, was interpolated by area =25 's mean stand-age value

South_AUS_empty_age <- subset(age_final, (area == 26 & new_age == 0) )
North_AUS_available_age <- subset(age_final, (area == 25 & new_age > 0))
South_AUS_empty_age$final_age <- mean(North_AUS_available_age$new_age)

all_age2 <- rbind(all_age,South_AUS_empty_age)

all_age3 <- all_age2[,c("lon","lat","final_age")]

#this points are filled by mean value
gg <- plot_map3(continent_age_fAPAR2, 
          varnam = "new_age",plot_title = "age before filled by AUS",
          latmin = -65, latmax = 85,combine=FALSE)

NPP_all <- read.csv("~/data/NPP_Yunke/NPP_Nmin_dataset_with_predictors.csv")
NPP_forest <- subset(NPP_all,pft=="Forest")

#show missing points and forest sites - this missing values have to be interpolated! Because otherwise, these forest plots will finally laed to NA points
gg$ggmap +
  geom_point(data=all_age3,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=NPP_forest,aes(lon,lat),col="blue",size=1.5)+
  theme_grey(base_size = 12)

#so our method here is to set mean values of each local continent, which is reasonable.


###below no need to run

#1. first check if this alternative method makes a difference - no

#if only merging AUS and don't deal with other plots
#5. merge new interpolated subset of age dataframe into last complete dataframe, and generating findal data to be used.
final <- merge(continent_age_fAPAR, South_AUS_empty_age[,c("lon","lat","final_age")], by=c("lon","lat"), all.x = T,sort = T)

final$age_used <- NA

final$age_used[is.na(final$final_age)==TRUE] <- final$new_age[is.na(final$final_age)==TRUE]
final$age_used[is.na(final$final_age)==FALSE] <- final$final_age[is.na(final$final_age)==FALSE]

final2 <- final[,c("lon","lat","age_used")]
names(final2) <- c("lon","lat","age")

final_data <- final2
final_data$age[final_data$age == 0] <- NA

summary(final_data)

#now, extracting site values from Tg, alpha, c/n.....
elev_nc <- read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc")
elev <- as.data.frame(nc_to_df(elev_nc, varnam = "elevation"))

final_age <- merge(final_data, elev, by=c("lon","lat"), all.x = T,sort = T)
names(final_age) <- c("lon","lat","age","z")

NPP_forest$new_mapped_age <- NA

library(spgwr)
a <- 1.5 # which degree (distance) of grid when interpolating gwr from global grids
#Extract Tg, PPFD, vpd, alpha,fAPAR,age,CNrt,LMA, max-vcmax25
for (i in 1:nrow(NPP_forest)) {
  tryCatch({
    #PPFD_total_fapar
    #age
    age_global <- na.omit(final_age)
    NRE_part <- subset(age_global,lon>(NPP_forest[i,"lon"]-a)&lon<(NPP_forest[i,"lon"]+a)&
                         lat>(NPP_forest[i,"lat"]-a)&lat<(NPP_forest[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- NPP_forest[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    NPP_forest[i,c("new_mapped_age")]  <- (gwr(age ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    print(i)
  }, error=function(e){})} 

NPP_forest$new_mapped_age[NPP_forest$new_mapped_age<=0] <- NA
summary(NPP_forest)

plot(NPP_forest$mapped_age~NPP_forest$new_mapped_age)

analyse_modobs2(NPP_forest,"mapped_age","new_mapped_age", type = "points",relative=TRUE)$gg 

aa <- subset(NPP_forest,is.na(new_mapped_age)==TRUE)

#show all map

vcmax25_df <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vcmax25.nc"),
  varnam = "vcmax25"))

Tg <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/Tg.nc"),
  varnam = "Tg"))

PPFD <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/PPFD.nc"),
  varnam = "PPFD"))

vpd <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vpd.nc"),
  varnam = "vpd"))

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

all_predictors_current <- as.data.frame(cbind(Tg$lon,Tg$lat,Tg$Tg, PPFD$PPFD, vpd$vpd,
                                      fAPAR$fAPAR, age$age,
                                      CNrt$CNrt, LMA$LMA, vcmax25_df$vcmax25 ))
dim(na.omit(all_predictors_current))

final_data2 <- final_data[order(final_data[,2],final_data[,1]),]

all_predictors_new <- as.data.frame(cbind(Tg$lon,Tg$lat,Tg$Tg, PPFD$PPFD, vpd$vpd,
                                              fAPAR$fAPAR, final_data2$age,
                                              CNrt$CNrt, LMA$LMA, vcmax25_df$vcmax25 ))
dim(na.omit(all_predictors_new))
55993-52105 # these will lose 3888 grids

#check where these "lost" plots located
all_predictors_old_new <- as.data.frame(cbind(Tg$lon,Tg$lat,Tg$Tg, PPFD$PPFD, vpd$vpd,
                                          fAPAR$fAPAR,age$age, final_data2$age,
                                          CNrt$CNrt, LMA$LMA, vcmax25_df$vcmax25 ))
all_predictors_old_new2 <- subset(all_predictors_old_new,is.na(V7)==F & is.na(V8)==T)
all_predictors_old_new3 <- all_predictors_old_new2[,c(1:7,9:11)]
all_predictors_old_new4 <- na.omit(all_predictors_old_new3)
dim(all_predictors_old_new4)
gg <- plot_map3(continent_age_fAPAR2, 
                varnam = "new_age",plot_title = "age before filled by AUS",
                latmin = -65, latmax = 85,combine=FALSE)

gg$ggmap +
  geom_point(data=all_predictors_old_new4,aes(V1,V2),col="red",size=1.5)+
  theme_grey(base_size = 12)
