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
devtools::load_all("/Users/yunpeng/yunkepeng/latest_packages/rbeni/") 
devtools::load_all("/Users/yunpeng/yunkepeng/latest_packages/ingestr/")
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)
library(spgwr)
##1.Sara Vicca
#1.1 NPP
Sara_NPP <- read.csv("~/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_NPP.csv")
Sara_NPP$Source_NPP <- Sara_NPP$Source
Sara_NPP <- Sara_NPP[,c("Plot","Begin.year","End.year","NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","ANPP_2","BNPP_1","TNPP_1","Source_NPP")]
#TNPP_1 = ANPP_2 + BNPP_1 (not including understory)
#ANPP_2 = NPP.foliage + NPP.wood

#convert NA begin year
Sara_NPP$Begin.year[Sara_NPP$Begin.year==9999] <- Sara_NPP$End.year[Sara_NPP$Begin.year==9999]
Sara_NPP$Begin.year[Sara_NPP$Begin.year==9999] <- Sara_NPP$End.year[Sara_NPP$Begin.year==9999]

#Sara_NPP$End.year - Sara_NPP$Begin.year

Sara_NPP$no <- c(1:nrow(Sara_NPP))

#1.2 siteinfo
Sara_NPP_siteinfo <- read.csv("~/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_siteinfo.csv")
Sara_NPP_siteinfo$Plot <- Sara_NPP_siteinfo$Plot.name

for (i in 1:nrow(Sara_NPP_siteinfo)){
  if (Sara_NPP_siteinfo$Direction.Latitude[i] == "S"){
    Sara_NPP_siteinfo$lat[i] <- -(Sara_NPP_siteinfo$Latitude[i])
  } else {
    Sara_NPP_siteinfo$lat[i] <- Sara_NPP_siteinfo$Latitude[i]
  }
  if (Sara_NPP_siteinfo$Direction.Longitude[i] == "W"){
    Sara_NPP_siteinfo$lon[i] <- -(Sara_NPP_siteinfo$Longitude[i])
  } else {
    Sara_NPP_siteinfo$lon[i] <- Sara_NPP_siteinfo$Longitude[i]
  }
}
Sara_NPP_siteinfo$Source_siteinfo <- Sara_NPP_siteinfo$Source.1
Sara_NPP_siteinfo <- Sara_NPP_siteinfo[,c("Plot","lon","lat","Elevation","Evergreen.Deciduous","Management.code","Management","Source_siteinfo")]

Sara_NPP2 <- merge(Sara_NPP,Sara_NPP_siteinfo,by=c("Plot"),all.x=TRUE)

#interpolate missing elevation - let's interpolate them by etopo
Sara_NPP2_elv_missing <- subset(Sara_NPP2,is.na(Elevation)==TRUE)
Sara_NPP2_elv_missing_Plot <- aggregate(Sara_NPP2_elv_missing,by=list(Sara_NPP2_elv_missing$Plot), FUN=mean, na.rm=TRUE)
Sara_NPP2_elv_missing_Plot <-Sara_NPP2_elv_missing_Plot[,c("Group.1","lon","lat")]
names(Sara_NPP2_elv_missing_Plot) <- c("sitename","lon","lat")
#devtools::load_all("~/yunkepeng/Grassland_new_ingestr_rsofun_20210326/ingestr/")
df_etopo <- ingest(
  Sara_NPP2_elv_missing_Plot,
  source = "etopo1",
  dir = "~/data/etopo/" 
)
Sara_NPP2_elv_missing_Plot$Elevation <- as.numeric(as.data.frame(df_etopo$data))
Sara_NPP2_elv_missing_Plot$Elevation
names(Sara_NPP2_elv_missing_Plot) <- c("Plot","lon","lat","Elevation_etopo")

#now, interpolate those NA elevation by etopo
Sara_NPP2 <- merge(Sara_NPP2,Sara_NPP2_elv_missing_Plot,by=c("Plot","lon","lat"),all.x=TRUE)
Sara_NPP2$Elevation[is.na(Sara_NPP2$Elevation)==TRUE] <- Sara_NPP2$Elevation_etopo[is.na(Sara_NPP2$Elevation)==TRUE]
Sara_NPP2$Elevation

#Now, add site-level stand-age
Sara_age <- read.csv("~/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_age.csv")
Sara_age <- Sara_age[,c("Plot","Stand.age")]
Sara_age_site <- aggregate(Sara_age,by=list(Sara_age$Plot), FUN=mean, na.rm=TRUE)
Sara_age_site <- subset(Sara_age_site,Stand.age>0)
Sara_age_site <- Sara_age_site[,c(1,3)]
names(Sara_age_site) <- c("Plot","age")
summary(Sara_age_site)

Sara_NPP3 <- merge(Sara_NPP2,Sara_age_site,by=c("Plot"),all.x=TRUE)

#Now, add site-level LAI
Sara_LAI <- read.csv("~/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_LAI.csv")
Sara_LAI <- Sara_LAI[,c("Plot","LAI")]
Sara_LAI_site <- aggregate(Sara_LAI,by=list(Sara_LAI$Plot), FUN=mean, na.rm=TRUE)
Sara_LAI_site <- Sara_LAI_site[,c(1,3)]
names(Sara_LAI_site) <- c("Plot","LAI")
Sara_LAI_site$observedfAPAR <- 1-exp(-0.5 * Sara_LAI_site$LAI)
hist(Sara_LAI_site$observedfAPAR)

Sara_NPP4 <- merge(Sara_NPP3,Sara_LAI_site,by=c("Plot"),all.x=TRUE)

#now, merged with GPP (1) primarily based on plot + start.year + end.year and (2) based on average of plot
#firstly, aggregate based on sitename, start.year and end.year
Sara_GPP <- read.csv("~/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_GPP.csv")
Sara_GPP <- Sara_GPP[,c("Plot","GPP","Begin.year","End.year")]
Sara_GPP$Begin.year[Sara_GPP$Begin.year==9999] <- Sara_GPP$End.year[Sara_GPP$Begin.year==9999]
Sara_GPP <- subset(Sara_GPP,GPP>0)
Sara_GPP_site <- aggregate(Sara_GPP,by=list(Sara_GPP$Plot,Sara_GPP$Begin.year,Sara_GPP$End.year), FUN=mean, na.rm=TRUE)
Sara_GPP_site <- Sara_GPP_site[,c(1,2,3,5)]
names(Sara_GPP_site) <- c("Plot","Begin.year","End.year","GPP")
Sara_NPP5 <- merge(Sara_NPP4,Sara_GPP_site,by=c("Plot","Begin.year","End.year"),all.x=TRUE)

#alternatively, aggregate based on site only
Sara_GPP <- read.csv("~/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_GPP.csv")
Sara_GPP <- Sara_GPP[,c("Plot","GPP")]
Sara_GPP_site <- aggregate(Sara_GPP,by=list(Sara_GPP$Plot), FUN=mean, na.rm=TRUE)
Sara_GPP_site <- subset(Sara_GPP_site,GPP>0)
Sara_GPP_site <- Sara_GPP_site[,c(1,3)]
names(Sara_GPP_site) <- c("Plot","GPP2")
Sara_NPP6 <- merge(Sara_NPP5,Sara_GPP_site,by=c("Plot"),all.x=TRUE)

for (i in 1:nrow(Sara_NPP6)){
  if (is.na(Sara_NPP6$GPP[i]) == TRUE){ 
    Sara_NPP6$GPP[i] <- Sara_NPP6$GPP2[i]
  } else {
    Sara_NPP6$GPP[i] <- Sara_NPP6$GPP[i]
  }
}

Sara_NPP6 <- Sara_NPP6[,!(names(Sara_NPP6) %in% "GPP2")]

#add alpha - as obtained earlier in SPLASH
#Sara_NPP6 <- Sara_NPP6[order(Sara_NPP6$no), ]
#alphalist3 <- read.csv(file="~/data/NPP_Yunke/NPP_Vicca/climates_alpha.csv")$alpha
#Sara_NPP6$alpha <- alphalist3

#add site-level soil C/N
Sara_CN <- read.csv(file="~/data/NPP_Yunke/NPP_Vicca/orig/References_Yunke_soilCN.csv")
Sara_CN <- Sara_CN[,c("Plot.name","Soil.C.N")]
Sara_CN$soilCN <- (as.numeric(gsub(",",".",Sara_CN[,2])))
hist(Sara_CN$soilCN)
Sara_CN_site <- aggregate(Sara_CN,by=list(Sara_CN$Plot.name), FUN=mean, na.rm=TRUE)
Sara_CN_site <- subset(Sara_CN_site,soilCN>0)
Sara_CN_site$Plot <- Sara_CN_site$Group.1
Sara_CN_site <- Sara_CN_site[,c("Plot","soilCN")]

Sara_NPP7 <- merge(Sara_NPP6,Sara_CN_site,by=c("Plot"),all.x=TRUE)
Sara_NPP7 <- Sara_NPP7[order(Sara_NPP7$no), ]
Sara_NPP7$pft <- "Forest"
Sara_NPP7$file <- "Sara Vicca"

#2. now, add Malhi # NPP_Malhi's data needs double check
NPP_Malhi <- read.csv("~/data/NPP_Yunke/NPP_Malhi/NPP_Malhi.csv")
NPP_Malhi <- NPP_Malhi[,c("site","lon","lat","z","file","Begin_year","End_year","Source","NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","ANPP_2","BNPP_1","TNPP_1","GPP")]
names(NPP_Malhi) <- c("Plot","lon","lat","Elevation","file","Begin.year","End.year","Source_NPP","NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","ANPP_2","BNPP_1","TNPP_1","GPP")
NPP_Malhi$Management.code <- "UM"
NPP_Malhi$pft <-"Forest"

#3. add Keith (take care about rep)
NPP_Keith <- read.csv("~/data/NPP_Yunke/NPP_Keith/orig/ABPE.csv")
NPP_Keith <- NPP_Keith[,c("Site","Ecosystem","age","lat","long","Elevation","Mgmt_code","ANPP","ANPP","GPP","Source")]
names(NPP_Keith) <- c("Plot","pft","age","lat","lon","Elevation","Management.code","ANPP_2","NPP.foliage","GPP","Source_NPP")
NPP_Keith$Begin.year <- 1991
NPP_Keith$End.year <- 2010
NPP_Keith$file <- "Keith"

NPP_Sara_Malhi_Keith <- dplyr::bind_rows(Sara_NPP7, NPP_Malhi,NPP_Keith) 

NPP_Sara_Malhi_Keith <- NPP_Sara_Malhi_Keith %>% 
  rename(
    site = Plot,
    Begin_year = Begin.year,
    End_year = End.year,
    z = Elevation)

head(NPP_Sara_Malhi_Keith)

##4. now, add ForC (take care about rep)
#1. Input original dataset of Forc
#NOT YET remove any FACE experiements (with co2 treatment)
#needs remove it at the end?

Forc <- read.csv(file="~/data/NPP_Yunke/NPP_ForC/orig/ForC_measurements.csv")
variablelist <- c("GPP_C","NPP_1_C","ANPP_woody_stem_C","ANPP_foliage_C","ANPP_woody_C","ANPP_2_C","BNPP_root_C","BNPP_root_coarse_C","BNPP_root_fine_C")
#according to definition in "~/data/NPP_Yunke/NPP_ForC/orig/ForC_variables.csv"
#NPP_1_C: NPP, including foliage, branch, stem, coarse root, and fine root. 
# ANPP_2_C (Annual aboveground NPP, including foliage, stem growth, and branch turnover)= ANPP_foliage_C + ANPP_woody_C:
#ANPP_woody_stem_C: stem ANPP, part of wood ANPP
# BNPP_root_C: Annual root production determined by repeated soil coring, isotopic estimates of fine-root turnover combined with biomass measurements, root ingrowth cores, root ingrowth mesh, upscaled root-length production observed in minirhizotrons or the soil respiration and litterfall constraint formulated by Raich & Nadelhoffer (1989)
length(variablelist)
mylist <- vector(mode = "list", length = length(variablelist))

for (i in 1:length(variablelist)){
  mylist[[i]] <- subset(Forc,variable.name==variablelist[i])
}

GPP_C <- as.data.frame(mylist[1])
NPP_1_C <- as.data.frame(mylist[2])
ANPP_woody_stem_C <- as.data.frame(mylist[3])
ANPP_foliage_C <- as.data.frame(mylist[4])
ANPP_woody_C <- as.data.frame(mylist[5])
ANPP_2_C <- as.data.frame(mylist[6])
BNPP_root_C <- as.data.frame(mylist[7])
BNPP_root_coarse_C <- as.data.frame(mylist[8])
BNPP_root_fine_C <- as.data.frame(mylist[9])

#firstly, find new GPP sites.
sites_NPP <- NPP_Sara_Malhi_Keith
sites_NPP$repeated <- 1
sites_NPP2 <- aggregate(repeated~site,data=sites_NPP,mean)
dim(sites_NPP2)

GPP_C <- as.data.frame(mylist[1])
GPP_C <- GPP_C[,c("sites.sitename","mean","date","start.date","end.date")]
names(GPP_C) <-c("site","GPP","date","start.date","end.date")

GNPP_final2 <- merge(GPP_C,sites_NPP2,by=c("site"),all.x=TRUE) # now, merge to get new GPP sites 
new_GPP <- subset(GNPP_final2,is.na(repeated)==TRUE)
new_GPP$year_start <- NA
new_GPP$year_end <- NA
for (i in 1:nrow(new_GPP)){
  if (is.na(new_GPP$start.date[i]) == FALSE){
    new_GPP$year_start[i] <- new_GPP$start.date[i]
    new_GPP$year_end[i] <- new_GPP$end.date[i]} else {
      new_GPP$year_start[i] <- new_GPP$date[i]
      new_GPP$year_end[i] <- new_GPP$date[i]}}

new_GPP$year_start <- round(as.numeric(new_GPP$year_start))  
new_GPP$year_end <- round(as.numeric(new_GPP$year_end))  
new_GPP$year_start[is.na(new_GPP$year_start)==TRUE] <- 1991
new_GPP$year_end[is.na(new_GPP$year_end)==TRUE] <- 2010
new_GPP <- new_GPP[,c("site","GPP","year_start","year_end")]
new_GPP$GPP<-100*new_GPP$GPP
hist(new_GPP$GPP)

#get npp dataset
#write a function to convert to correct measurement yaer - and aggrgated based on mean~sites.sitename+year_start+year_end. So that it can be merged later on
object_correct_years <- function(object_npp){
  object_npp$year_start <- NA
  object_npp$year_end <- NA
  for (i in 1:nrow(object_npp)){
    if (is.na(object_npp$start.date[i]) == FALSE){
      object_npp$year_start[i] <- object_npp$start.date[i]
      object_npp$year_end[i] <- object_npp$end.date[i]} else {
        object_npp$year_start[i] <- object_npp$date[i]
        object_npp$year_end[i] <- object_npp$date[i]}}
  
  object_npp$year_start <- round(as.numeric(object_npp$year_start))  
  object_npp$year_end <- round(as.numeric(object_npp$year_end))  
  object_npp$year_start[is.na(object_npp$year_start)==TRUE] <- 1991
  object_npp$year_end[is.na(object_npp$year_end)==TRUE] <- 2010
  object_npp <- object_npp[,c("sites.sitename","year_start","year_end","mean")]
  object_npp$mean <- object_npp$mean * 100 #convert to consistent unit gC/m2/yr
  object_npp2 <- aggregate(mean~sites.sitename+year_start+year_end,data=object_npp,mean,na.rm=TRUE)
  return(object_npp2)
}

#with additional info - including ambient/elevated info
object_correct_years_with_siteinfo <- function(object_npp){
  object_npp$year_start <- NA
  object_npp$year_end <- NA
  for (i in 1:nrow(object_npp)){
    if (is.na(object_npp$start.date[i]) == FALSE){
      object_npp$year_start[i] <- object_npp$start.date[i]
      object_npp$year_end[i] <- object_npp$end.date[i]} else {
        object_npp$year_start[i] <- object_npp$date[i]
        object_npp$year_end[i] <- object_npp$date[i]}}
  
  object_npp$year_start <- round(as.numeric(object_npp$year_start))  
  object_npp$year_end <- round(as.numeric(object_npp$year_end))  
  object_npp$year_start[is.na(object_npp$year_start)==TRUE] <- 1991
  object_npp$year_end[is.na(object_npp$year_end)==TRUE] <- 2010
  object_npp <- object_npp[,c("sites.sitename","plot.name","year_start","year_end","mean")]
  object_npp$mean <- object_npp$mean * 100 #convert to consistent unit gC/m2/yr
  object_npp2 <- aggregate(mean~sites.sitename+plot.name+year_start+year_end,data=object_npp,mean,na.rm=TRUE)
  return(object_npp2)
}

NPP_1_C <- object_correct_years_with_siteinfo(NPP_1_C)
ANPP_woody_stem_C <- object_correct_years(ANPP_woody_stem_C)
ANPP_foliage_C<- object_correct_years(ANPP_foliage_C)
ANPP_woody_C<- object_correct_years(ANPP_woody_C)
ANPP_2_C<- object_correct_years(ANPP_2_C)
BNPP_root_C<- object_correct_years(BNPP_root_C)
BNPP_root_coarse_C<- object_correct_years(BNPP_root_coarse_C)
BNPP_root_fine_C<- object_correct_years(BNPP_root_fine_C)

names(NPP_1_C) <- c("sites.sitename","Management","year_start","year_end","TNPP_1")
names(ANPP_woody_stem_C) <- c("sites.sitename","year_start","year_end","NPP.stem")
names(ANPP_foliage_C) <- c("sites.sitename","year_start","year_end","NPP.foliage")
names(ANPP_woody_C) <- c("sites.sitename","year_start","year_end","NPP.wood")
names(ANPP_2_C) <- c("sites.sitename","year_start","year_end","ANPP_2")
names(BNPP_root_C) <- c("sites.sitename","year_start","year_end","BNPP_1")
names(BNPP_root_coarse_C) <- c("sites.sitename","year_start","year_end","NPP.coarse")
names(BNPP_root_fine_C) <- c("sites.sitename","year_start","year_end","NPP.fine")


NPP_all <-Reduce(function(x,y) merge(x = x, y = y, by = c("sites.sitename","year_start","year_end"),all.x=TRUE),
                 list(NPP_1_C,ANPP_woody_stem_C,ANPP_foliage_C,ANPP_woody_C,ANPP_2_C,BNPP_root_C,BNPP_root_coarse_C,BNPP_root_fine_C))
names(sites_NPP2) <- c("sites.sitename","repeated")
NPP_final2 <- merge(NPP_all,sites_NPP2,by=c("sites.sitename"),all.x=TRUE) # now, merge to get new NPP sites 
new_NPP <- subset(NPP_final2,is.na(repeated)==TRUE)
dim(new_NPP)

#now, merged new_GPP to new_NPP - primarily based on sitename+year
names(new_GPP) <- c("sites.sitename","GPP_sitename_yr","year_start","year_end")
new_GPP_site_yr <- aggregate(GPP_sitename_yr~sites.sitename+year_start+year_end,new_GPP,FUN=mean, na.rm=TRUE)

new_NPP_GPP1 <- merge(new_NPP,new_GPP_site_yr,by=c("sites.sitename","year_start","year_end"),all.x=TRUE) # now, merge to get new GPP sites 
summary(new_NPP_GPP1) #merged with 23 points - ok

#And secondly based on sitename only.
names(new_GPP) <- c("sites.sitename","GPP_sitename","year_start","year_end")
new_GPP_sitemean <- aggregate(GPP_sitename~sites.sitename,new_GPP,FUN=mean, na.rm=TRUE)
new_NPP_GPP2 <- merge(new_NPP_GPP1,new_GPP_sitemean,by=c("sites.sitename"),all.x=TRUE) # now, merge to get new GPP sites 

#primarily using GPP_sitename_yr, and secondily using GPP_sitename
new_NPP_GPP3 <- new_NPP_GPP2 %>% mutate(GPP = coalesce(GPP_sitename_yr,GPP_sitename))

#now, removing new_GPP sites
ForC_all <- new_NPP_GPP3[,!(names(new_NPP_GPP3) %in% c("repeated","GPP_sitename_yr","GPP_sitename"))]

#finally, combine with coordinates
forc_coord <- read.csv("~/data/NPP_Yunke/NPP_ForC/orig/ForC_sites.csv")
forc_coord <- forc_coord[,c("sites.sitename","lat","lon","masl")]
ForC_all_coord <- merge(ForC_all,forc_coord,by=c("sites.sitename"),all.x=TRUE) # now, merge to get new GPP sites 
summary(ForC_all_coord)
# a few sites have missed coordinates - needs combined manually, based on ForC_sites.csv above
ForC_all_coord$lat[ForC_all_coord$sites.sitename=="Cascade Head 1 "] <- 45.1024
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="Cascade Head 1 "]<- -123.8816
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="Cascade Head 1 "]<- 205

ForC_all_coord$lat[ForC_all_coord$sites.sitename=="Deer Canyon Preserve Pinyon Juniper Woodland "]<- 34.36
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="Deer Canyon Preserve Pinyon Juniper Woodland "]<--106.27
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="Deer Canyon Preserve Pinyon Juniper Woodland "]<-2126

ForC_all_coord$lat[ForC_all_coord$sites.sitename=="O site "] <- 44.5
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="O site "] <- -121.617
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="O site "] <- 915

ForC_all_coord$lat[ForC_all_coord$sites.sitename=="Tablelands Juniper Savanna "]<-34.43
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="Tablelands Juniper Savanna "]<- -105.86
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="Tablelands Juniper Savanna "] <- 1926

ForC_all_coord$lat[ForC_all_coord$sites.sitename=="University of Michigan Biological Station (UMBS) "]<- 45.583
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="University of Michigan Biological Station (UMBS) "]<- -84.7
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="University of Michigan Biological Station (UMBS) "]<- NA

ForC_all_coord$lat[ForC_all_coord$sites.sitename=="Valles Caldera Mixed Conifer "] <- 35.89
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="Valles Caldera Mixed Conifer "] <- -106.53
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="Valles Caldera Mixed Conifer "] <- 3049

ForC_all_coord$lat[ForC_all_coord$sites.sitename=="Valles Caldera Ponderosa Pine "] <- 35.86
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="Valles Caldera Ponderosa Pine "] <- -106.6
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="Valles Caldera Ponderosa Pine "] <- 2486

ForC_all_coord$lat[ForC_all_coord$sites.sitename=="Willow Creek (WC)-Chequamegon National Forest "]<- 45.783
ForC_all_coord$lon[ForC_all_coord$sites.sitename=="Willow Creek (WC)-Chequamegon National Forest "]<- -90.083
ForC_all_coord$masl[ForC_all_coord$sites.sitename=="Willow Creek (WC)-Chequamegon National Forest "]<- 480

ForC_all_coord$masl <- as.numeric(ForC_all_coord$masl)

ForC_elv_missing <- subset(ForC_all_coord,is.na(masl)==TRUE)
ForC_elv_missing_Plot <- aggregate(ForC_elv_missing,by=list(ForC_elv_missing$sites.sitename), FUN=mean, na.rm=TRUE)
ForC_elv_missing_Plot <-ForC_elv_missing_Plot[,c("Group.1","lon","lat")]
names(ForC_elv_missing_Plot) <- c("sitename","lon","lat")
head(ForC_elv_missing_Plot)

df_etopo <- ingest(
  ForC_elv_missing_Plot,
  source = "etopo1",
  dir = "~/data/etopo/" 
)
ForC_elv_missing_Plot$Elevation <- as.numeric(as.data.frame(df_etopo$data))
names(ForC_elv_missing_Plot) <- c("sites.sitename","lon","lat","Elevation_etopo")

#now, interpolate those NA elevation by etopo
ForC_all_coord2 <- merge(ForC_all_coord,ForC_elv_missing_Plot,by=c("sites.sitename","lon","lat"),all.x=TRUE)
ForC_all_coord2$Elevation <- ForC_all_coord2$masl
ForC_all_coord2$Elevation[is.na(ForC_all_coord2$masl)==TRUE] <- ForC_all_coord2$Elevation_etopo[is.na(ForC_all_coord2$masl)==TRUE]

ForC_all_coord2 <- ForC_all_coord2[,!(names(ForC_all_coord2) %in% c("masl","Elevation_etopo"))]
ForC_all_coord2 <- ForC_all_coord2 %>% 
  rename(site=sites.sitename,
         Begin_year=year_start,
         End_year=year_end,
         z=Elevation)

#aggregate basing lon, lat, z, begin_year, end_year (this is important because some data is flipped from the same thing - but it doesn't change)
#ForC_all_coord2_final <- aggregate(ForC_all_coord2,by=list(ForC_all_coord2$lon,ForC_all_coord2$lat,ForC_all_coord2$z,ForC_all_coord2$Begin_year,ForC_all_coord2$End_year,ForC_all_coord2$site,ForC_all_coord2$Management), mean,na.rm=TRUE)
#ForC_all_coord2_final$site <- ForC_all_coord2_final$Group.6
#ForC_all_coord2_final$Management <- ForC_all_coord2_final$Group.7
#ForC_all_coord2_final <- ForC_all_coord2_final[,!(names(ForC_all_coord2_final) %in% c("Group.1","Group.2","Group.3","Group.4","Group.5","Group.6","Group.7"))]
#dim(ForC_all_coord2_final)

ForC_all_coord2$pft <- "Forest"
ForC_all_coord2$file <- "ForC"
ForC_all_coord2$Source_NPP <- "ForC"

NPP_Sara_Malhi_Keith_Forc <- dplyr::bind_rows(NPP_Sara_Malhi_Keith, ForC_all_coord2) 

NPP_all <- NPP_Sara_Malhi_Keith_Forc[,!(names(NPP_Sara_Malhi_Keith_Forc) %in% c("no","Source_siteinfo","Elevation_etopo"))]

summary(NPP_all)

#5. add Schulze 
NPP_Schulze <- read.csv(file="~/data/NPP_Yunke/NPP_Schulze/NPP_Schulze.csv")
NPP_Schulze <- NPP_Schulze[,c("site","lon","lat","z","file","Begin_year","End_year","NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","ANPP_2","BNPP_1","TNPP_1")]
NPP_Schulze$Management.code <- "UM"
NPP_Schulze$pft <-"Forest"
NPP_Schulze$file <- "NPP_Schulze"
NPP_Schulze$Source_NPP <- "NPP_Schulze"

NPP_all <- dplyr::bind_rows(NPP_all, NPP_Schulze) 


## 6. add data from Tian Di (pft = grassland for all data)
#firstly, clean our current data
Tiandi_df <- read.csv(file="~/data/npp_stoichiometry_grasslands_tiandi/npp_stoichiometry_china_grassland_CN_stoichiometry_with_matched_NPP_data_from_Prof_Fang_group_20201026.csv")
#as proved in Beni's ref, there is no big diff about lon_stoichmenistry and lon_npp, so we used the lon_npp because it is npp analyses now!
Tiandi_npp <- Tiandi_df[,c("Original_Site_Label_stoichiometry","Longitude_stoichiometry","Latitude_stoichiometry","Altitude_stoichiometry","Sample_time_NPP","Sample_time_NPP","TNPP","ANPP","ANPP","BNPP","CNratio_leaf","CNratio_root","CNratio_soil")]
#we will go back to c/n ratio later! It is NPP now only...
names(Tiandi_npp) <- c("site","lon","lat","z","Begin_year","End_year","TNPP_1","ANPP_2","NPP.foliage","BNPP_1","CN_leaf","CN_root","CN_soil")

#correct measurement year!
for (i in 1:nrow(Tiandi_npp)){
  if (is.na(Tiandi_npp$Begin_year[i]) == TRUE){ #if measruement year not available
    Tiandi_npp$Begin_year[i] <- 1991#convert to long-term
    Tiandi_npp$End_year[i] <- 2010 #convert to long-term 
  } else {
    Tiandi_npp$Begin_year[i] <- as.numeric(substr(Tiandi_npp$Begin_year[i], start = 1, stop = 4))
    Tiandi_npp$End_year[i] <- as.numeric(substr(Tiandi_npp$End_year[i], start = nchar(Tiandi_npp$End_year[i])-3, stop = nchar(Tiandi_npp$End_year[i]))) #1-4 or 6-9
  }
}

Tiandi_npp$Begin_year <- as.numeric(Tiandi_npp$Begin_year)
Tiandi_npp$End_year <- as.numeric(Tiandi_npp$End_year)
summary(Tiandi_npp$Begin_year)
summary(Tiandi_npp$End_year)
summary(Tiandi_npp)
Tiandi_npp$file <- "Tiandi Grassland"
Tiandi_npp$pft <- "Grassland"
Tiandi_npp$Source_NPP <- "Tiandi Grassland"

## 7. add data from Campioli (pft = grassland for all data)
Cam_df <- read.csv(file="~/data/campioli/grasslands_MCampioli_20160111.csv")
#correct coordinates firstly
for (i in 1:nrow(Cam_df)){
  if (Cam_df$latitude_sign[i] == "S"){
    Cam_df$lat[i] <- -(Cam_df$latitude_value[i])
  } else {
    Cam_df$lat[i] <- Cam_df$latitude_value[i]
  }
  if (Cam_df$longitude_sign[i] == "W"){
    Cam_df$lon[i] <- -(Cam_df$longitude_value[i])
  } else {
    Cam_df$lon[i] <- Cam_df$longitude_value[i]
  }
}

Cam_npp <- Cam_df[,c("site","lon","lat","elevation","period_start","period_end","tnpp","anpp","anpp","bnpp","managment","biome")]
names(Cam_npp) <- c("site","lon","lat","z","Begin_year","End_year","TNPP_1","ANPP_2","NPP.foliage","BNPP_1","management_MCampioli","biome_MCampioli")

#rbind them and manually add some input
Cam_npp$file <- "MCampioli"
Cam_npp$pft <- Cam_npp$biome_MCampioli
Cam_npp$Management.code <- Cam_npp$management_MCampioli
Cam_npp$Source_NPP <- "MCampioli"

NPP_final <- dplyr::bind_rows(NPP_all, Tiandi_npp,Cam_npp) 

summary(NPP_final)
dim(NPP_final)

NPP_final <- NPP_final[,!(names(NPP_final) %in% c("management_MCampioli","biome_MCampioli","CN_soil"))]

## 8. add more forest sites from corrected Sara Vicca's dataset, including anpp, npp.leaf and npp.wood. (leafcn data will be included below)
Sara2_df <- read.csv(file="~/data/NPP_Yunke/NPP_Vicca/orig/CORRECTIONS_CascadeHead_Andrews.csv")
Sara2_df2 <- subset(Sara2_df,Repeat=="no") #remove repeated data as inputted in NPP_SaraVicca

Sara2_NPP <- Sara2_df2[,c("LONGITUDE","LATITUDE","ELEVATION","YEAR","YEAR","AG_PROD_TREE_TOTAL_AS_CARBON","AG_PROD_TREE_FOLIAGE_AS_CARBON","AG_PROD_TREE_WOOD_AS_CARBON")]
#unit of AG_PROD_TREE_TOTAL_AS_CARBON and others: gC/m2/yr
names(Sara2_NPP) <- c("lon","lat","z","Begin_year","End_year","ANPP_2","NPP.foliage","NPP.wood")
#remove values == -9999 and 0, defined as NA
Sara2_NPP$ANPP_2[Sara2_NPP$ANPP_2<=0] <- NA
Sara2_NPP$NPP.foliage[Sara2_NPP$NPP.foliage<=0] <- NA
Sara2_NPP$NPP.wood[Sara2_NPP$NPP.wood <=0] <- NA
Sara2_NPP$file <- "Vicca_validation_file"
Sara2_NPP$pft<-"Forest"
Sara2_NPP$Source_NPP <- "Sara Vicca"

summary(Sara2_NPP)

#create site name
Sara2_NPP_sitename <- aggregate(Sara2_NPP,by=list(Sara2_NPP$lon,Sara2_NPP$lat,Sara2_NPP$z), mean,na.rm=TRUE)
for (i in 1:nrow(Sara2_NPP_sitename)){
  Sara2_NPP_sitename$site[i] <- paste("Sara2_NPP",i,sep = "")
}

Sara2_NPP_sitename <- Sara2_NPP_sitename[,c("lon","lat","z","site")]

Sara2_NPP$no <- c(1:nrow(Sara2_NPP))

Sara2_NPP2 <- merge(Sara2_NPP,Sara2_NPP_sitename,by=c("lon","lat","z"),all.x=TRUE)

Sara2_NPP2 <- Sara2_NPP2[order(Sara2_NPP2$no), ]

Sara2_NPP2 <- Sara2_NPP2[,!(names(Sara2_NPP2) %in% "no")]

NPP_final2 <- dplyr::bind_rows(NPP_final, Sara2_NPP2) 
summary(NPP_final2)

###Now, merge leaf c/n from various sources

# (1) Add Schulz - unit: all in g/m2

CN_Schulz <- read.csv(file="~/data/NPP_Yunke/npp_cn/CN_Schulze.csv")
CN_Schulz2 <- CN_Schulz[,c(5,48:58)]

CN_Schulz2$CN_leaf_Schulz <- CN_Schulz2$c_leaf/CN_Schulz2$n_leaf
CN_Schulz2$CN_root_Schulz <- (CN_Schulz2$c_coarseroot+CN_Schulz2$c_fineroot)/CN_Schulz2$n_root
#???(CN_Schulz2$c_fineroot + CN_Schulz2$c_coarseroot) / CN_Schulz2$n_root

CN_Schulz2$CN_stem_Schulz <- CN_Schulz2$c_stem/CN_Schulz2$n_stem
CN_Schulz2$CN_wood_Schulz <- (CN_Schulz2$c_stem+CN_Schulz2$c_branch)/(CN_Schulz2$n_stem+CN_Schulz2$n_branch)

CN_Schulz2 <- CN_Schulz2[,c("site","CN_leaf_Schulz","CN_root_Schulz","CN_stem_Schulz","CN_wood_Schulz")]
CN_Schulz2

# (2) Add Malhi data 
CN_Malhi <- read.csv(file="~/data/NPP_Yunke/npp_cn/CN_Malhi.csv")

CN_Malhi$CN_leaf_alt_malhi <- CN_Malhi$cmass/CN_Malhi$nmass 
CN_Malhi <- subset(CN_Malhi,CN_leaf_alt_malhi>0)
CN_Malhi2 <- CN_Malhi[,c("site","CN_leaf_alt_malhi")]

#The step is (see comment below):
#(1) aggregate leaf C/N based on lon + lat, according to a leaf traits dataset provided by Sara Vicca. There are more than one 1 record in a site, beacuse it measured many trees (individuals) in one site. See their original csv.
CN_SaraVicca <- read.csv(file="~/data/NPP_Yunke/NPP_Vicca/orig/NACP_TERRA_PNW_leaf_trait.csv")
CN_SaraVicca <- CN_SaraVicca[,c("LONGITUDE","LATITUDE","LEAF_CN")]
CN_SaraVicca$LEAF_CN[CN_SaraVicca$LEAF_CN<0] <- NA
CN_SaraVicca$LONGITUDE <- as.numeric(CN_SaraVicca$LONGITUDE)
CN_SaraVicca$LATITUDE <- as.numeric(CN_SaraVicca$LATITUDE)

CN_SaraVicca2 <- aggregate(CN_SaraVicca, by=list(CN_SaraVicca$LONGITUDE,CN_SaraVicca$LATITUDE), mean,na.rm=TRUE)
CN_SaraVicca2 <- CN_SaraVicca2[,c("LONGITUDE","LATITUDE","LEAF_CN")]
names(CN_SaraVicca2) <- c("lon","lat","CN_leaf_sara")
dim(CN_SaraVicca2)
summary(CN_SaraVicca2)

#(2) Merge C/N to NPP dataset STRICTLY based on lon and lat, which nearly all merged to org file with (CORRECT2..)
test <- merge(NPP_final2,CN_SaraVicca2,by=c("lon","lat"),all.x=TRUE)
test1 <- subset(test,CN_leaf_sara>0)
test1 <- test1[,c("lon","lat","z","CN_leaf_sara")]
dim(test1) # make them aggregate basing on lon + lat + z, to prepare for merge then
test1 <- aggregate(test1, by=list(test1$lon,test1$lat,test1$z), mean,na.rm=TRUE)
dim(test1)
test1<- test1[,c("lon","lat","z","CN_leaf_sara")]

#(3) now, merge with NPP_final2 
#firstly merge based on site, for CN_Schulze and CN_Malhi
npp_cn1 <-Reduce(function(x,y) merge(x = x, y = y, by = c("site"),all.x=TRUE),
                 list(NPP_final2,CN_Schulz2,CN_Malhi2))

#merge based on test1 (New Sara's dataset - for strictly merged condition)
npp_cn2 <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat","z"),all.x=TRUE),
                 list(npp_cn1,test1))

npp_cn2$CN_stem_final <- npp_cn2$CN_stem_Schulz
npp_cn2$CN_wood_final <- npp_cn2$CN_wood_Schulz

#now, it is the time to combine
npp_cn3 <- npp_cn2 %>% mutate(CN_leaf_final = coalesce(CN_leaf,CN_leaf_Schulz,CN_leaf_alt_malhi,CN_leaf_sara)) %>%
  mutate(CN_root_final = coalesce(CN_root,CN_root_Schulz))


hist(npp_cn3$CN_leaf_final)
hist(npp_cn3$CN_stem_final)
hist(npp_cn3$CN_wood_final)
hist(npp_cn3$CN_root_final)
NPP_final3 <- npp_cn3[,!(names(npp_cn3) %in% c("CN_leaf","CN_root","CN_leaf_Schulz","CN_root_Schulz","CN_stem_Schulz","CN_wood_Schulz","CN_leaf_alt_malhi","CN_leaf_sara"))]

NPP_final3$lnf_obs_final <-NPP_final3$NPP.foliage/NPP_final3$CN_leaf_final
NPP_final3$bnf_obs_final  <- NPP_final3$BNPP_1/NPP_final3$CN_root_final
NPP_final3$wnf_obs_final  <- NPP_final3$NPP.wood/NPP_final3$CN_wood_final

summary(NPP_final3) 

## 9. add Tiandi's latest forest data (the last data we need!)
tiandi_forest <- read.csv(file="~/data/npp_stoichiometry_forests_tiandi/Site_level_forest_CN_NPP_China_TD_20210104_for_Beni_Yunke.csv")
tiandi_forest <- tiandi_forest[,c(1,2,3,4,11,17,18,19,20,21,22)]
head(tiandi_forest)
names(tiandi_forest) <- c("site","lon","lat","z","CN_leaf_final","CN_root_final","TNPP_1","ANPP_2",
                          "NPP.foliage","NPP.wood","BNPP_1")

#convert unit from tC/ha/yr to gC/m2/yr --> *100
tiandi_forest$TNPP_1 <- 100*tiandi_forest$TNPP_1
tiandi_forest$ANPP_2 <- 100*tiandi_forest$ANPP_2
tiandi_forest$NPP.foliage <- 100*tiandi_forest$NPP.foliage
tiandi_forest$NPP.wood <- 100*tiandi_forest$NPP.wood
tiandi_forest$BNPP_1 <- 100*tiandi_forest$BNPP_1
tiandi_forest$file <- "~/data/npp_stoichiometry_forests_tiandi/"
tiandi_forest$Begin_year <- 2006
tiandi_forest$End_year <- 2015
tiandi_forest$pft <- "Forest"

tiandi_forest$lnf_obs_final <-tiandi_forest$NPP.foliage/tiandi_forest$CN_leaf_final 
tiandi_forest$bnf_obs_final  <- tiandi_forest$BNPP_1/tiandi_forest$CN_root_final

NPP_final7 <- dplyr::bind_rows(NPP_final3, tiandi_forest) 

#one last thing - correct a few Malhi's site lon/lat
#Correct Malhi's wrong point (their paper's original lon/lat were mixed! so we correct them here - by turning them around)
NPP_final7$lon[NPP_final7$site=="BCI Plateau, Panama?"] <- -79.85
NPP_final7$lat[NPP_final7$site=="BCI Plateau, Panama?"] <- 9.17

NPP_final7$lon[NPP_final7$site=="BDFFP Fazenda"] <- -60.17
NPP_final7$lat[NPP_final7$site=="BDFFP Fazenda"] <-  -2.63

NPP_final7$lon[NPP_final7$site=="Bionte, Brazil"] <- -60.17
NPP_final7$lat[NPP_final7$site=="Bionte, Brazil"] <- -2.63

NPP_final7$lon[NPP_final7$site=="Mocambo, Brazil"] <- -48.45
NPP_final7$lat[NPP_final7$site=="Mocambo, Brazil"] <- -1.45

NPP_final7$lon[NPP_final7$site== "San Carlos caatinga"] <- -67.05
NPP_final7$lat[NPP_final7$site== "San Carlos caatinga"] <- 1.75

NPP_final7$lon[NPP_final7$site== "San Carlos terra firme"] <- -67.05
NPP_final7$lat[NPP_final7$site== "San Carlos terra firme"] <- 1.93

NPP_final7$lon[NPP_final7$site=="Tapajo?s, Brazil" & NPP_final7$z==100] <- -55
NPP_final7$lat[NPP_final7$site=="Tapajo?s, Brazil"& NPP_final7$z==100] <- -2.75

#convert pft
NPP_final7$pft[NPP_final7$pft=="grassland"] <- "Grassland"

NPP <- NPP_final7

NPP %>% group_by(pft)  %>% summarise(number = n())
NPP$NPP.foliage[NPP$NPP.foliage==0] <- NA
NPP$NPP.stem[NPP$NPP.stem==0] <- NA
NPP$NPP.wood[NPP$NPP.wood==0] <- NA
NPP$NPP.fine[NPP$NPP.fine==0] <- NA
NPP$NPP.coarse[NPP$NPP.coarse==0] <- NA
NPP$ANPP_2[NPP$ANPP_2==0] <- NA
NPP$TNPP_1[NPP$TNPP_1==0] <- NA


#remove repeated data
NPP$rep <- "not_repeated"

#remove rep data 
NPP$rep[NPP$site=="Waring's Woods"] <- "repeated" #npp, bnpp, anpp all repeated

#remove rep data from Keith, where repeated to Campioli or Vicca, the definition of repeatation is:
#below are those samples with the same site-name to Campioli, Vicca or ForC dataset
#because Keith's dataset was actually derived from Sara Vicca - we don't want repeat here - if Keith and others file has the same site-name (or coordinates data), remove data from Keith and only keep others data
NPP$rep[NPP$site=="CA-Let-F01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="CG-tch-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="CN-Inn-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="CN-Inn-F01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="DE-gri-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="KZ-shr-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="RU-ha1-F01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="RU-ha2-F01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="RU-ha3-F01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="RU-krs-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="US-jas-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="US-kbs-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="US-kon-D05"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="US-osg-D01"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Morgan Monroe"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Bornhoved Alder"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Bornhoved Beech"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Davos"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Flagstaff, USA"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Tapajos, Brazil"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Oregon, USA"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Metolius"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Caxiuana, Brazil"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Coastal plain North Carolina"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Turkey Point TP74"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Turkey Point TP89"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Turkey Point TP02"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="UCI-1998"&NPP$file=="Keith"] <- "repeated"
NPP$rep[NPP$site=="Takayama"&NPP$file=="Keith"] <- "repeated"

#finally, all forest plots in Keith was removed, only keep a few grassland sites

#remove rep data (paired repeated measurements between ForC and Sara Vicca) 
NPP$rep[NPP$site=="Thompson dry chronosequence d12"&NPP$file=="ForC"] <- "repeated"
NPP$rep[NPP$site=="Thompson dry chronosequence d131"&NPP$file=="ForC"] <- "repeated"
NPP$rep[NPP$site=="Thompson dry chronosequence d20"&NPP$file=="ForC"] <- "repeated"
NPP$rep[NPP$site=="Thompson dry chronosequence d37"&NPP$file=="ForC"] <- "repeated"
NPP$rep[NPP$site=="Thompson dry chronosequence d71"&NPP$file=="ForC"] <- "repeated"
NPP$rep[NPP$site=="Willow Creek (WC)-Chequamegon National Forest"&NPP$file=="ForC"] <- "repeated"
NPP$rep[NPP$site=="Wind River Canopy Crane"&NPP$file=="ForC"] <- "repeated"
NPP$rep[NPP$site=="Qianyanzhou"&NPP$file=="ForC"] <- "repeated"

#removing part for collecting c3,c4 (too confusing and not used in this submission). For info see submission/Forest_site_org.R 

#additional step: not yet confirmed - so not included it now - and gpp will not be used in this dataset
#add more gpp data from Compioli et al. SI table 1 
#NPP$GPP[NPP$site=="IT-bea-D02"] <- 1568
#NPP$GPP[NPP$site=="US-che-D01"] <- 626
#NPP$GPP[NPP$site=="DE-gri-D01"] <- 1233
#NPP$GPP[NPP$site=="CN-Hab-F01"] <- 634
#NPP$GPP[NPP$site=="RU-ha1-F01"] <- 519
#NPP$GPP[NPP$site=="RU-ha3-F01"] <- 526
#NPP$GPP[NPP$site=="CN-Inn-D01_C"] <- 182
#NPP$GPP[NPP$site=="US-kbs-D01"] <- 1015
#NPP$GPP[NPP$site=="US-kbs-D04"] <- 512
#NPP$GPP[NPP$site=="US-kbs-D05"] <- 374
#NPP$GPP[NPP$site=="US-kbs-D03"] <- 793
#NPP$GPP[NPP$site=="US-jas-D01"] <- 516
#NPP$GPP[NPP$site=="US-kon-D05"] <- 1151
#NPP$GPP[NPP$site=="RU-krs-D01"] <- 1611
#NPP$GPP[NPP$site=="CA-Let-F01"] <- 280
#NPP$GPP[NPP$site=="CA-mat-D01"] <- 786
#NPP$GPP[NPP$site=="US-osg-D01"] <- 1890
#NPP$GPP[NPP$site=="CG-tch-D01"] <- 1572
#NPP$GPP[NPP$site=="US-Spe-D01"] <- 829
#add more gpp data from keith's source (with same site name)
#NPP$GPP[NPP$site=="CN-Inn-F01"] <- 204 

##### Finally Input N uptake
#(1) newly added Nmin rate data from Finzi
Finzi <- read.csv("~/data/NPP_Yunke/Nmin_Finzi/Nmin_Finzi.csv")
names(Finzi)[names(Finzi) == "Lat"] <- "lat"
names(Finzi)[names(Finzi) == "Long"] <- "lon"

#Forest - only merging forest this time
Finzi_Forest <- subset(Finzi, Biome!="temp grass")
unique(Finzi_Forest$Biome)
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

Nmin_dataset <- Nmin_final[,c("lon","lat","elv","sitename","year_start","year_end","Publication","Nmin")]
names(Nmin_dataset) <- c("lon","lat","z","site","Begin_year","End_year","Source_NPP","Nmin")
Nmin_dataset$pft<- "Forest"
Nmin_dataset$file<- "Finzi"
Nmin_dataset$rep<- "not_repeated"

NPP_Nuptake <- dplyr::bind_rows(NPP, Nmin_dataset) 

NPP_Nuptake$year_start <- NPP_Nuptake$Begin_year
NPP_Nuptake$year_end <- NPP_Nuptake$End_year

NPP_Nuptake$year_start[NPP_Nuptake$Begin_year<=1980] <- 1980
NPP_Nuptake$year_end[NPP_Nuptake$End_year<=1980] <- 1989

csvfile <- paste("~/data/NPP_Yunke/NPP_Nmin_dataset.csv")
write_csv(NPP_Nuptake, path = csvfile)

#please note! below 3 steps requires above output NPP_Nmin_dataset.csv
#but I guess, This csv in the future may be changed/updated, if some filter more required
#However, as far as lon, lat, z, begin_year and end_year from (NPP_Nmin_dataset.csv) not changed/added. It's fine to not replicate the (too long) process below. Just directly merging with below output

#combine with p-model's gpcessing/pmodel_simulation.R
gpp_vcmax25 <- read.csv("~/data/NPP_Yunke/simulated_gpp/site_simulated_gpp_vcmax.csv")
NPP_Nuptake_gpp_vcmax25 <- merge(NPP_Nuptake,gpp_vcmax25,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)

#combine site-simulated vpd, ppfd, alpha and Tg: see climate_site_data.R
climates_sites <- read.csv("~/data/NPP_Yunke/predictors/climates_sites.csv")
NPP_Nuptake_gpp_vcmax25_climates <- merge(NPP_Nuptake_gpp_vcmax25,climates_sites,by=c("lon","lat","z","year_start","year_end","Begin_year","End_year"),all.x=TRUE)
summary(NPP_Nuptake_gpp_vcmax25_climates)

#combine with site-interpolated value from all maps: see map_site_data.R
gwr_sites <- read.csv("~/data/NPP_Yunke/predictors/predictors_gwr.csv")
gwr_sites <- gwr_sites %>% 
  rename(
    mapped_age = age)

NPP_Nuptake_gpp_vcmax25_climates_gwr <- merge(NPP_Nuptake_gpp_vcmax25_climates,gwr_sites,by=c("lon","lat","z","year_start","year_end","Begin_year","End_year"),all.x=TRUE)
summary(NPP_Nuptake_gpp_vcmax25_climates_gwr)


#final dataset - now conducts consistency and repeated chceck
dataset1 <- subset(NPP_Nuptake_gpp_vcmax25_climates_gwr,rep=="not_repeated") #removed repeated

dataset1[grep("elevated", dataset1$Management),]$rep <- "elevated_co2" #remove with co2 elevated measurements
dataset1[grep("Elevated", dataset1$Management),]$rep <- "elevated_co2" 

#remove co2 data!
dataset2 <- subset(dataset1,rep=="not_repeated")

#convert grassland npp.foliage = anpp
#dataset2$NPP.foliage[dataset2$pft!="Forest"] <- dataset2$ANPP_2[dataset2$pft!="Forest"]

#check npp = anpp + bnpp
dataset2$check_npp <- dataset2$ANPP_2+dataset2$BNPP_1-dataset2$TNPP_1
aa <-subset(dataset2,check_npp< -5 |check_npp >5 )[c(1:8,14:16,27,52)] 
dim(aa)
dataset2$outlier[dataset2$check_npp< -5 |dataset2$check_npp >5  ] <- "NPP!=ANPP+BNPP"
#check anpp = npp.foliage + npp.wood
dataset2$check_anpp <- dataset2$NPP.foliage+dataset2$NPP.wood-dataset2$ANPP_2
bb <-subset(dataset2,check_anpp< -5 |check_anpp>5 )[c(1:8,9,11,14,27,53)]
dim(bb)
dataset2$outlier[dataset2$check_anpp< -5 |dataset2$check_anpp >5  ] <- "ANPP!=NPP.foliage+NPP.wood"

#they need figured out
###Make sure anpp= wood +leaf; npp = anpp + bnpp
dataset2$ANPP_2[is.na(dataset2$NPP.foliage)==FALSE&is.na(dataset2$NPP.wood)==FALSE] <- dataset2$NPP.foliage[is.na(dataset2$NPP.foliage)==FALSE&is.na(dataset2$NPP.wood)==FALSE] +dataset2$NPP.wood[is.na(dataset2$NPP.foliage)==FALSE&is.na(dataset2$NPP.wood)==FALSE]
dataset2$TNPP_1[is.na(dataset2$ANPP_2)==FALSE&is.na(dataset2$BNPP_1)==FALSE] <- dataset2$ANPP_2[is.na(dataset2$ANPP_2)==FALSE&is.na(dataset2$BNPP_1)==FALSE] +dataset2$BNPP_1[is.na(dataset2$ANPP_2)==FALSE&is.na(dataset2$BNPP_1)==FALSE]

summary(dataset2$ANPP_2+dataset2$BNPP_1-dataset2$TNPP_1)
summary(dataset2$NPP.foliage+dataset2$NPP.wood-dataset2$ANPP_2)

#check if it is balanced and without special NA
# npp.foliage is missing, but has anpp and wood - filling them
subset(dataset2,is.na(ANPP_2)==F &is.na(NPP.wood)==F &pft=="Forest"&is.na(NPP.foliage)==T)
dataset2$NPP.foliage[is.na(dataset2$ANPP_2)==F &is.na(dataset2$NPP.wood)==F &dataset2$pft=="Forest"&is.na(dataset2$NPP.foliage)==T] <- 
  dataset2$ANPP_2[is.na(dataset2$ANPP_2)==F &is.na(dataset2$NPP.wood)==F &dataset2$pft=="Forest"&is.na(dataset2$NPP.foliage)==T] -
  dataset2$NPP.wood[is.na(dataset2$ANPP_2)==F &is.na(dataset2$NPP.wood)==F &dataset2$pft=="Forest"&is.na(dataset2$NPP.foliage)==T]

# npp.wood is missing, but has anpp and npp.foliage - filling them
subset(dataset2,is.na(ANPP_2)==F &is.na(NPP.foliage)==F &pft=="Forest"&is.na(NPP.wood)==T)
dataset2$NPP.wood[is.na(dataset2$ANPP_2)==F &is.na(dataset2$NPP.foliage)==F &dataset2$pft=="Forest"&is.na(dataset2$NPP.wood)==T] <- 
  dataset2$ANPP_2[is.na(dataset2$ANPP_2)==F &is.na(dataset2$NPP.foliage)==F &dataset2$pft=="Forest"&is.na(dataset2$NPP.wood)==T] -
  dataset2$NPP.foliage[is.na(dataset2$ANPP_2)==F &is.na(dataset2$NPP.foliage)==F &dataset2$pft=="Forest"&is.na(dataset2$NPP.wood)==T]

# anpp is missing, but has NPP and BNPP - filling them
subset(dataset2,is.na(TNPP_1)==F &is.na(BNPP_1)==F &is.na(ANPP_2)==T)
dataset2$ANPP_2[is.na(dataset2$TNPP_1)==F &is.na(dataset2$BNPP_1)==F &is.na(dataset2$ANPP_2)==T] <- 
  dataset2$TNPP_1[is.na(dataset2$TNPP_1)==F &is.na(dataset2$BNPP_1)==F &is.na(dataset2$ANPP_2)==T] -
  dataset2$BNPP_1[is.na(dataset2$TNPP_1)==F &is.na(dataset2$BNPP_1)==F &is.na(dataset2$ANPP_2)==T]

# bnpp is missing, but has NPP and ANPP - filling them
subset(dataset2,is.na(TNPP_1)==F &is.na(ANPP_2)==F &is.na(BNPP_1)==T)
dataset2$BNPP_1[is.na(dataset2$TNPP_1)==F &is.na(dataset2$ANPP_2)==F &is.na(dataset2$BNPP_1)==T] <- 
  dataset2$TNPP_1[is.na(dataset2$TNPP_1)==F &is.na(dataset2$ANPP_2)==F &is.na(dataset2$BNPP_1)==T] -
  dataset2$ANPP_2[is.na(dataset2$TNPP_1)==F &is.na(dataset2$ANPP_2)==F &is.na(dataset2$BNPP_1)==T]

#double check
summary(dataset2$ANPP_2+dataset2$BNPP_1-dataset2$TNPP_1)
summary(dataset2$NPP.foliage+dataset2$NPP.wood-dataset2$ANPP_2)

#have observed that some repeated data ForC. In this way, remove such repeated data
dataset2$addtional_rep <- duplicated(dataset2[,c("lon","lat","z","Begin_year","End_year","TNPP_1","ANPP_2","BNPP_1")])
dim(subset(dataset2,addtional_rep=="TRUE" & file=="ForC"))

dataset2$rep_ForC[dataset2$addtional_rep=="TRUE" & dataset2$file=="ForC"] <- "repeated"
dataset3 <- subset(dataset2,is.na(rep_ForC)==TRUE)
dim(dataset2)
dim(dataset3)

dataset4 <- dataset3[,!(names(dataset3) %in% c("check_npp","check_anpp","outlier","check_bnpp","addtional_rep","rep_ForC","rep"))]

#only keep forest and grassland
dataset5 <- subset(dataset4,pft=="Forest" | pft=="Grassland")

#one thing
dataset5$NPP.foliage[dataset5$site=="Sara2_NPP72"]
dataset5$ANPP_2[dataset5$site=="Sara2_NPP72"]
dataset5$NPP.foliage[dataset5$site=="Sara2_NPP72"]  <- NA #one thing shown wrong npp.foliage, with npp.leaf = anpp, remove npp.leaf

#remove the old file from NPP_Schulze
dataset6 <- subset(dataset5,file!="NPP_Schulze")

#check leaf and root C/N
summary(aggregate(subset(dataset6,pft=="Grassland"),by=list(subset(dataset6,pft=="Grassland")$site), FUN=mean, na.rm=TRUE)$CN_leaf_final)
summary(aggregate(subset(dataset6,pft=="Grassland"),by=list(subset(dataset6,pft=="Grassland")$site), FUN=mean, na.rm=TRUE)$CN_root_final)

#check how many sites (in lon and lat or site - no big difference! just use site-defined version)
dim(unique(subset(dataset6,pft=="Grassland" & CN_leaf_final>0)[,c("lon","lat")]))
length(unique(subset(dataset6,pft=="Grassland" & CN_leaf_final>0)$site))

dim(unique(subset(dataset6,pft=="Grassland" & CN_root_final>0)[,c("lon","lat")]))
length(unique(subset(dataset6,pft=="Grassland" & CN_root_final>0)$site))

#REMOVE Tiandi Grassland's npp and bnpp
dataset6 <- subset(dataset6,file!="Tiandi Grassland")

#convert grassland stand-age as NA (e.g. Keith: from 1 to NA)
dataset6$age[dataset6$pft=="Grassland"] <- NA

#final calculation of lnf, bnf and wnf 
dataset6$lnf_obs_final <-dataset6$NPP.foliage/dataset6$CN_leaf_final
dataset6$bnf_obs_final  <- dataset6$BNPP_1/dataset6$CN_root_final
dataset6$wnf_obs_final  <- dataset6$NPP.wood/dataset6$CN_wood_final


#remove columns not used
dataset6 <- dataset6[,!(names(dataset6) %in% c("Evergreen.Deciduous","Management.code",
                                               "Management","pred_gpp_c3","max_vcmax25_c3",
                                               "alpha_sites","PPFD_sites","Tg_sites","vpd_sites",
                                               "PPFD_total_fapar","PPFD_total","alpha","mapped_gpp"))]
#check if all stem < wood ->yes 
summary(dataset6$NPP.stem/dataset6$NPP.wood)
subset(dataset6,NPP.stem+NPP.foliage-ANPP_2>0) # this three samples from "ORNL-FACE" are strange - where stem+foliage - ANPP_2 >0 - we removed stem values here since it is not reasonable
dataset6$NPP.stem[dataset6$site=="ORNL-FACE"] <- NA
  
#check npp.fine or npp.coarse
summary(dataset6$NPP.fine+dataset6$NPP.coarse-dataset6$BNPP_1)
unique(subset(dataset6,NPP.fine+NPP.coarse-BNPP_1>5 |NPP.fine+NPP.coarse-BNPP_1< -5)$site) # these site's corase and fine root converted to NA
dataset6$NPP.fine[dataset6$site=="Podocarpus 1"] <- NA;dataset6$NPP.coarse[dataset6$site=="Podocarpus 1"] <- NA
dataset6$NPP.fine[dataset6$site=="ORNL-FACE"] <- NA;dataset6$NPP.coarse[dataset6$site=="ORNL-FACE"] <- NA
dataset6$NPP.fine[dataset6$site=="POP-EUROFACE"] <- NA;dataset6$NPP.coarse[dataset6$site=="POP-EUROFACE"] <- NA
dataset6$NPP.fine[dataset6$site=="Xiaohu NF"] <- NA;dataset6$NPP.coarse[dataset6$site=="Xiaohu NF"] <- NA
dataset6$NPP.fine[dataset6$site=="Jadraas Class I"] <- NA;dataset6$NPP.coarse[dataset6$site=="Jadraas Class I"] <- NA
summary(dataset6$NPP.fine+dataset6$NPP.coarse-dataset6$BNPP_1) # else are tiny - it is fine

csvfile <- paste("~/data/NPP_Yunke/NPP_Nmin_dataset_with_predictors.csv")
write_csv(dataset6, path = csvfile)

#check: plot missing data - p model's vcmax25 - many of them are missing due to on the edge

#already checked all missing sites, comparing with past attempt of validations in NPP_final or NPP_grassland_final, all consistent.

#this seems be due to fapar in those places are NA, try re-working on this plots in rsofun to fill vcmax25?
aa <- subset(NPP_Nuptake_gpp_vcmax25_climates_gwr,is.na(max_vcmax25_c3)==TRUE)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(aa$lon,aa$lat, col="red", pch=16,cex=1)

#check: plot missing data - mapped vcmax25 - many of them are missing due to on the edge
aa <- subset(NPP_Nuptake_gpp_vcmax25_climates_gwr,is.na(vcmax25)==TRUE)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(aa$lon,aa$lat, col="red", pch=16,cex=1)


#part 2, prepare NRE dataset
#2. Input NRE from different sources, rbind, remove repeated data and plot maps, and test a lm model
#Dong is species-based, Du is site-based. To make them consistent, let's use site-based for both df

NRE_Du <- read.csv(file="~/data/NRE_various/NRE_Du/NRE_Du.csv")
NRE_Dong <- read.csv(file="~/data/NRE_various/NRE_Deng/NRE_Deng.csv")

#first - make forest model only
NRE_Du_df <- NRE_Du[,c("lon","lat","NRE","MAT","MAP","VegeType")]
#vegetype =1 is woody ecosystem (all assumed as forest here)
#vegetype = 2 is grassland ecosystem
NRE_Du_df <- subset(NRE_Du_df,VegeType==1)

NRE_Du_df <- aggregate(NRE_Du_df,by=list(NRE_Du_df$lon,NRE_Du_df$lat), FUN=mean, na.rm=TRUE) #site-mean
NRE_Du_df <- NRE_Du_df[,c(3:7)]
head(NRE_Du_df)
dim(NRE_Du_df)

NRE_Dong_df <- NRE_Dong[,c("Longitude","Latitude","NRE.nitrogen.resorption.efficiency.","MAT","MAP","Biome.abbreviation...")]
names(NRE_Dong_df) <- c("lon","lat","NRE","MAT","MAP","biome")
NRE_Dong_df <- subset(NRE_Dong_df,biome=="TRF"|biome=="STF"|biome=="TF"|biome=="BF")
#Forest is TRF,STF,TF,BF
#Desert is Des; Tundra is TUN
#Grassland is Grs
NRE_Dong_df <- aggregate(NRE_Dong_df,by=list(NRE_Dong_df$lon,NRE_Dong_df$lat), FUN=mean, na.rm=TRUE) #site-mean
head(NRE_Dong_df)
NRE_Dong_df <- NRE_Dong_df[,c(3:7)]
dim(NRE_Dong_df)


NRE_Dong_df$source <- "Dong"
NRE_Du_df$source <- "Du"
NRE_df <- rbind(NRE_Du_df,NRE_Dong_df)
summary(NRE_df)

#check repeated data, and remove 4 repeated points from Du et al. paper
NRE_df$repeated <- duplicated(NRE_df[,c("lon","lat")])
subset(NRE_df,repeated==TRUE)
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

#4. Extract site climate/soil/age data (from prediction fields in nimpl simulation) for all NRE sites
#Input data from: ~/data/nimpl_sofun_inputs/map/Final_ncfile
library(rbeni)

#input elevation for global grids
elev_nc <- read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc")
elev <- as.data.frame(nc_to_df(elev_nc, varnam = "elevation"))
head(elev) 

#input nc file
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

#cbind all predictors, and its lon, lat, z
all_predictors <- cbind(elev,Tg$Tg,PPFD$PPFD,vpd$vpd,
                        alpha$alpha,fAPAR$fAPAR)

names(all_predictors) <- c("lon","lat","z","Tg","PPFD","vpd",
                           "alpha","fAPAR")

Tg_df <- all_predictors[,c("lon","lat","z","Tg")]
PPFD_df <- all_predictors[,c("lon","lat","z","PPFD")]
vpd_df <- all_predictors[,c("lon","lat","z","vpd")]
alpha_df <- all_predictors[,c("lon","lat","z","alpha")]
fAPAR_df <- all_predictors[,c("lon","lat","z","fAPAR")]

#now, apply gwr to extract site predictors' value
head(NRE_df)
NRE_df$elevation[is.na(NRE_df$elevation)==TRUE] <- 0

NRE_site <- NRE_df[,c("lon","lat","elevation")]
names(NRE_site) <- c("lon","lat","z")


a <- 1.5 # which degree (distance) of grid when interpolating gwr from global grids
i <- 1
#Extract Tg, PPFD, vpd, alpha
for (i in c(1:nrow(NRE_site))){ #one site does not have elevation (NA), therefore omitted
  #Tg
  Tg_global <- na.omit(Tg_df)
  NRE_part <- subset(Tg_global,lon>(NRE_site[i,1]-a)&lon<(NRE_site[i,1]+a)&
                       lat>(NRE_site[i,2]-a)&lat<(NRE_site[i,2]+a))
  coordinates(NRE_part) <- c("lon","lat")
  gridded(NRE_part) <- TRUE
  NRE_coord <- NRE_site[i,1:3]
  coordinates(NRE_coord) <- c("lon","lat")
  NRE_site$Tg[i] <- (gwr(Tg ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
  #ppfd
  PPFD_global <- na.omit(PPFD_df)
  NRE_part <- subset(PPFD_global,lon>(NRE_site[i,1]-a)&lon<(NRE_site[i,1]+a)&
                       lat>(NRE_site[i,2]-a)&lat<(NRE_site[i,2]+a))
  coordinates(NRE_part) <- c("lon","lat")
  gridded(NRE_part) <- TRUE
  NRE_coord <- NRE_site[i,1:3]
  coordinates(NRE_coord) <- c("lon","lat")
  NRE_site$PPFD[i] <- (gwr(PPFD ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
  #vpd
  vpd_global <- na.omit(vpd_df)
  NRE_part <- subset(vpd_global,lon>(NRE_site[i,1]-a)&lon<(NRE_site[i,1]+a)&
                       lat>(NRE_site[i,2]-a)&lat<(NRE_site[i,2]+a))
  coordinates(NRE_part) <- c("lon","lat")
  gridded(NRE_part) <- TRUE
  NRE_coord <- NRE_site[i,1:3]
  coordinates(NRE_coord) <- c("lon","lat")
  NRE_site$vpd[i] <- (gwr(vpd ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
  #alpha
  alpha_global <- na.omit(alpha_df)
  NRE_part <- subset(alpha_global,lon>(NRE_site[i,1]-a)&lon<(NRE_site[i,1]+a)&
                       lat>(NRE_site[i,2]-a)&lat<(NRE_site[i,2]+a))
  coordinates(NRE_part) <- c("lon","lat")
  gridded(NRE_part) <- TRUE
  NRE_coord <- NRE_site[i,1:3]
  coordinates(NRE_coord) <- c("lon","lat")
  NRE_site$alpha[i] <- (gwr(alpha ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
}

#NRE_site newly including Ndep
library(hwsdr)
#devtools::load_all("/Users/yunpeng/yunkepeng/compuetational_ingestr/ingestr/")

NRE_site$nhx <- NA
NRE_site$noy <- NA

#reference of Du and Deng mostly from 1991-2009
for (i in 1:nrow(NRE_site)) {
  tryCatch({
    print(i)
    df_ndep <- ingest_bysite(
      sitename  = paste("a",i,sep=""),
      source    = "ndep",
      lon       = NRE_site$lon[i],
      lat       = NRE_site$lat[i],
      year_start= 1991,
      year_end  = 2009,
      timescale = "y",
      dir       = "~/data/ndep_lamarque/",
      verbose   = FALSE
    )
    NRE_site$noy[i] <- mean(df_ndep$noy,na.rm=TRUE)
    NRE_site$nhx[i] <- mean(df_ndep$nhx,na.rm=TRUE)
  }, error=function(e){})} 

NRE_site$ndep <- NRE_site$noy + NRE_site$nhx

#combine NRE with site extracted climate and CNrt (Tg, PPFD, vpd, alpha, CNrt)
NRE_climate <- cbind(NRE_df[,c("NRE","MAT","MAP","source")],NRE_site)
summary(NRE_climate)

NRE_climate$vpd[NRE_climate$vpd<=0] <-NA
NRE_climate$nre <- NRE_climate$NRE/100
summary(NRE_climate)
csvfile <- paste("~/data/NRE_various/NRE_dataset.csv")
write.csv(NRE_climate, csvfile, row.names = TRUE)
