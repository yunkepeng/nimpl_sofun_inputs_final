rm(list=ls)
##1.Sara Vicca
#1.1 NPP
Sara_NPP <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_NPP.csv")
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
Sara_NPP_siteinfo <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_siteinfo.csv")
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
devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/ingestr/")
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
Sara_age <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_age.csv")
Sara_age <- Sara_age[,c("Plot","Stand.age")]
Sara_age_site <- aggregate(Sara_age,by=list(Sara_age$Plot), FUN=mean, na.rm=TRUE)
Sara_age_site <- subset(Sara_age_site,Stand.age>0)
Sara_age_site <- Sara_age_site[,c(1,3)]
names(Sara_age_site) <- c("Plot","age")
summary(Sara_age_site)

Sara_NPP3 <- merge(Sara_NPP2,Sara_age_site,by=c("Plot"),all.x=TRUE)

#Now, add site-level LAI
Sara_LAI <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_LAI.csv")
Sara_LAI <- Sara_LAI[,c("Plot","LAI")]
Sara_LAI_site <- aggregate(Sara_LAI,by=list(Sara_LAI$Plot), FUN=mean, na.rm=TRUE)
Sara_LAI_site <- Sara_LAI_site[,c(1,3)]
names(Sara_LAI_site) <- c("Plot","LAI")
Sara_LAI_site$observedfAPAR <- 1-exp(-0.5 * Sara_LAI_site$LAI)
hist(Sara_LAI_site$observedfAPAR)

Sara_NPP4 <- merge(Sara_NPP3,Sara_LAI_site,by=c("Plot"),all.x=TRUE)

#now, merged with GPP (1) primarily based on plot + start.year + end.year and (2) based on average of plot
#firstly, aggregate based on sitename, start.year and end.year
Sara_GPP <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_GPP.csv")
Sara_GPP <- Sara_GPP[,c("Plot","GPP","Begin.year","End.year")]
Sara_GPP$Begin.year[Sara_GPP$Begin.year==9999] <- Sara_GPP$End.year[Sara_GPP$Begin.year==9999]
Sara_GPP <- subset(Sara_GPP,GPP>0)
Sara_GPP_site <- aggregate(Sara_GPP,by=list(Sara_GPP$Plot,Sara_GPP$Begin.year,Sara_GPP$End.year), FUN=mean, na.rm=TRUE)
Sara_GPP_site <- Sara_GPP_site[,c(1,2,3,5)]
names(Sara_GPP_site) <- c("Plot","Begin.year","End.year","GPP")
Sara_NPP5 <- merge(Sara_NPP4,Sara_GPP_site,by=c("Plot","Begin.year","End.year"),all.x=TRUE)

#alternatively, aggregate based on site only
Sara_GPP <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/Forests_Colin_GPP.csv")
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
Sara_NPP6 <- Sara_NPP6[order(Sara_NPP6$no), ]
alphalist3 <- read.csv(file="/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/climates_alpha.csv")$alpha
Sara_NPP6$alpha <- alphalist3

#add site-level soil C/N
Sara_CN <- read.csv(file="/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/References_Yunke_soilCN.csv")
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

#statistical model
library(lme4)
library(nlme)
library(lmerTest)
library("PerformanceAnalytics")
library(MuMIn)
library(tidyverse)
mod_anpp <- lmer(log((ANPP_2/GPP)/(1-(ANPP_2/GPP))) ~ log(soilCN) + log(age) + alpha + observedfAPAR + (1|Plot), data = Sara_NPP7)
summary(mod_anpp)
r.squaredGLMM(mod_anpp)

mod_tnpp <- lmer(log((TNPP_1/GPP)/(1-(TNPP_1/GPP))) ~ log(soilCN) + log(age) + observedfAPAR  + (1|Plot), data = Sara_NPP7)
summary(mod_tnpp)
r.squaredGLMM(mod_tnpp)

#2. now, add Malhi
NPP_Malhi <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Malhi/NPP_Malhi.csv")
NPP_Malhi <- NPP_Malhi[,c("site","lon","lat","z","file","Begin_year","End_year","Source","NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","ANPP_2","BNPP_1","TNPP_1","GPP")]
names(NPP_Malhi) <- c("Plot","lon","lat","Elevation","file","Begin.year","End.year","Source_NPP","NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","ANPP_2","BNPP_1","TNPP_1","GPP")
NPP_Malhi$Management.code <- "UM"
NPP_Malhi$pft <-"Forest"

#3. add Keith (take care about rep)
NPP_Keith <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Keith/orig/ABPE.csv")
NPP_Keith <- NPP_Keith[,c("Site","Ecosystem","age","lat","long","Elevation","Mgmt_code","ANPP","GPP","Source")]
names(NPP_Keith) <- c("Plot","pft","age","lat","lon","Elevation","Management.code","ANPP_2","GPP","Source_NPP")
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

Forc <- read.csv(file="~/data/NPP_Yunke/NPP_ForC/orig/ForC_measurements.csv")
variablelist <- c("GPP_C","NPP_1_C","ANPP_woody_stem_C","ANPP_foliage_C","ANPP_woody_C","ANPP_2_C","BNPP_root_C","BNPP_root_coarse_C","BNPP_root_fine_C")
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

NPP_1_C <- object_correct_years(NPP_1_C)
ANPP_woody_stem_C <- object_correct_years(ANPP_woody_stem_C)
ANPP_foliage_C<- object_correct_years(ANPP_foliage_C)
ANPP_woody_C<- object_correct_years(ANPP_woody_C)
ANPP_2_C<- object_correct_years(ANPP_2_C)
BNPP_root_C<- object_correct_years(BNPP_root_C)
BNPP_root_coarse_C<- object_correct_years(BNPP_root_coarse_C)
BNPP_root_fine_C<- object_correct_years(BNPP_root_fine_C)

names(NPP_1_C) <- c("sites.sitename","year_start","year_end","TNPP_1")
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
NPP_final2 <- merge(NPP_all,sites_NPP2,by=c("sites.sitename"),all.x=TRUE) # now, merge to get new GPP sites 
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
forc_coord <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_ForC/orig/ForC_sites.csv")
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
names(ForC_all_coord2) <- c("site","lon","lat","Begin_year","End_year","TNPP_1","NPP.stem","NPP.foliage","NPP.wood","ANPP_2","BNPP_1","NPP.coarse","NPP.fine","GPP","z")
ForC_all_coord2$pft <- "Forest"
ForC_all_coord2$file <- "ForC"

NPP_Sara_Malhi_Keith_Forc <- dplyr::bind_rows(NPP_Sara_Malhi_Keith, ForC_all_coord2) 

summary(NPP_Sara_Malhi_Keith_Forc)

NPP_all <- NPP_Sara_Malhi_Keith_Forc[,!(names(NPP_Sara_Malhi_Keith_Forc) %in% c("no","Source_siteinfo","Elevation_etopo"))]

#5. add Schulze 
NPP_Schulze <- read.csv(file="~/data/NPP_Yunke/NPP_Schulze/NPP_Schulze.csv")
NPP_Schulze <- NPP_Schulze[,c("site","lon","lat","z","file","Begin_year","End_year","Source","NPP.foliage","NPP.stem","NPP.wood","NPP.fine","NPP.coarse","ANPP_2","BNPP_1","TNPP_1")]
NPP_Schulze$Management.code <- "UM"
NPP_Schulze$pft <-"Forest"
NPP_Schulze$file <- "NPP_Schulze"

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

# XXX Beni:
Tiandi_npp %>% 
  ggplot(aes(x = CN_leaf, y = ..count..)) +
  geom_histogram()

#XXX: These values are much lower than the ones in the Terra-P dataset. We discussed that this may be a real difference. To be sure, could you please check with Di Tian whether the data here is also in units of gC / gN?
#YYY: Yes, exactly. But I have already checked with Tian Di that the data is gC/gN. Not sure why it is smaller, will investigate their reasons further. 

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
Tiandi_npp$file <- "Tiandi Grassland"
Tiandi_npp$pft <- "Grassland"

Cam_npp$file <- "MCampioli"
Cam_npp$pft <- Cam_npp$biome_MCampioli
Cam_npp$Management.code <- Cam_npp$management_MCampioli

NPP_final <- dplyr::bind_rows(NPP_all, Tiandi_npp,Cam_npp) 

summary(NPP_final)
dim(NPP_final)

NPP_final <- NPP_final[,!(names(NPP_final) %in% c("management_MCampioli","biome_MCampioli","CN_soil"))]

## 8. add more forest sites from corrected Sara Vicca's dataset, including anpp, npp.leaf and npp.wood. (leafcn data will be included below)
Sara2_df <- read.csv(file="~/data/NPP_Yunke/NPP_SaraVicca/orig/validation_data/CORRECTIONS_CascadeHead_Andrews.csv")
Sara2_df2 <- subset(Sara2_df,Repeat=="no") #remove repeated data as inputted in NPP_SaraVicca
Sara2_NPP <- Sara2_df2[,c("LONGITUDE","LATITUDE","ELEVATION","YEAR","YEAR","AG_PROD_TREE_TOTAL_AS_CARBON","AG_PROD_TREE_FOLIAGE_AS_CARBON","AG_PROD_TREE_WOOD_AS_CARBON")]
names(Sara2_NPP) <- c("lon","lat","z","Begin_year","End_year","ANPP_2","NPP.foliage","NPP.wood")
#remove values == -9999 and 0, defined as NA
Sara2_NPP$ANPP_2[Sara2_NPP$ANPP_2<=0] <- NA
Sara2_NPP$NPP.foliage[Sara2_NPP$NPP.foliage<=0] <- NA
Sara2_NPP$NPP.wood[Sara2_NPP$NPP.wood <=0] <- NA
Sara2_NPP$file <- "~/data/NPP_Yunke/NPP_SaraVicca/orig/validation_data"
Sara2_NPP$pft<-"Forest"
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

#XXX is this m2 ground area? Please specify units in the README.
#YYY Done

#XXX Distinction between stem and branches appears quite important (very different C:N ratios!). Is this disinction made in other datasets? What does ‘wood’ represent in the other datasets?
#Yes, agree that stem and branch is quite important here with different C/N so I have included them both. For wood or stem C/N ratio, this is the ONLY dataset we hold.
CN_Schulz <- read.csv(file="~/data/NPP_Yunke/npp_cn/CN_Schulze.csv")
CN_Schulz2 <- CN_Schulz[,c(5,48:58)]

CN_Schulz2$CN_leaf_Schulz <- CN_Schulz2$c_leaf/CN_Schulz2$n_leaf
CN_Schulz2$CN_root_Schulz <- (CN_Schulz2$c_coarseroot+CN_Schulz2$c_fineroot)/CN_Schulz2$n_root
#???(CN_Schulz2$c_fineroot + CN_Schulz2$c_coarseroot) / CN_Schulz2$n_root

CN_Schulz2$CN_stem_Schulz <- CN_Schulz2$c_stem/CN_Schulz2$n_stem
CN_Schulz2$CN_wood_Schulz <- (CN_Schulz2$c_stem+CN_Schulz2$c_branch)/(CN_Schulz2$n_stem+CN_Schulz2$n_branch)

CN_Schulz2 <- CN_Schulz2[,c("site","CN_leaf_Schulz","CN_root_Schulz","CN_stem_Schulz","CN_wood_Schulz")]
CN_Schulz2

# (2) Add Malhi data - assume cmass as constant 0.48 g/g; narea in gm-2, lma in gm-2
CN_Malhi <- read.csv(file="~/data/NPP_Yunke/npp_cn/CN_Malhi.csv")

#No original data of cmass but we can assume cmass = 46%, because (1) it is consistent with what we find in mean values of ~/data/leaf_traits/combined_leaf_traits.csv, equals to 46% and (2) see Enquist et al. 2017 https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12645 - fig.2, overall the cmass was within a very small variance through this elevation transect, and we can just assume this value as 0.48!
#set Average Cmass as 46%, based on latest (and largest) leaf traits database in ~/data/leaf_traits/combined_leaf_traits.csv

CN_Malhi$CN_leaf_alt_malhi <- CN_Malhi$cmass/CN_Malhi$nmass 
CN_Malhi <- subset(CN_Malhi,CN_leaf_alt_malhi>0)
CN_Malhi2 <- CN_Malhi[,c("site","CN_leaf_alt_malhi")]


#(3) Add Species-based traits data, as provided from Sara Vicca, including around 40 forest sites that have species-based leaf c/n
#The data is c% and n%

#XXX add variable descriptions for file NACP_TERRA_PNW_leaf_trait.csv in the README, sitting in the same directory.
#YYY: Done

# YYY: The step is (see comment below):
#(1) aggregate leaf C/N based on lon + lat, according to a leaf traits dataset provided by Sara Vicca. There are more than one 1 record in a site, beacuse it measured many trees (individuals) in one site. See their original csv.
CN_SaraVicca <- read.csv(file="~/data/NPP_Yunke/NPP_SaraVicca/orig/validation_data/NACP_TERRA_PNW_leaf_trait.csv")
CN_SaraVicca <- CN_SaraVicca[,c("LONGITUDE","LATITUDE","LEAF_CN")]
CN_SaraVicca$LEAF_CN[CN_SaraVicca$LEAF_CN<0] <- NA
CN_SaraVicca$LONGITUDE <- as.numeric(CN_SaraVicca$LONGITUDE)
CN_SaraVicca$LATITUDE <- as.numeric(CN_SaraVicca$LATITUDE)

CN_SaraVicca2 <- aggregate(CN_SaraVicca, by=list(CN_SaraVicca$LONGITUDE,CN_SaraVicca$LATITUDE), mean,na.rm=TRUE)
CN_SaraVicca2 <- CN_SaraVicca2[,c("LONGITUDE","LATITUDE","LEAF_CN")]
names(CN_SaraVicca2) <- c("lon","lat","CN_leaf_sara")
dim(CN_SaraVicca2)
summary(CN_SaraVicca2)

#(2) Merge C/N to NPP dataset STRICTLY based on lon and lat
test <- merge(NPP_final2,CN_SaraVicca2,by=c("lon","lat"),all.x=TRUE)
dim(subset(test,CN_leaf_sara>0))
test1 <- subset(test,CN_leaf_sara>0)
test1 <- test1[,c("lon","lat","z","CN_leaf_sara")]
dim(test1)
test1 <- aggregate(test1, by=list(test1$lon,test1$lat,test1$z), mean,na.rm=TRUE)
dim(test1)
test1<- test1[,c("lon","lat","z","CN_leaf_sara")]

#(3) Merge C/N to NPP dataset based on lon and lat, but not strictly
#--> we assume that it can be best merged within 0.01 lon/lat resolution
# I would persist using this way because it could help to generate more leaf NPP * C/N sites (from Sara and ForC), it would help us generate more validation sites 
NPP_old <- subset(NPP_final,file=="ForC"|file=="Sara Vicca") #only to select "ForC" and "Sarra Vicca" to be merged with more leaf C/N
NPP_old$lon <- round(NPP_old$lon,2)
NPP_old$lat <- round(NPP_old$lat,2)

CN_SaraVicca3 <- CN_SaraVicca
names(CN_SaraVicca3) <- c("lon","lat","CN_SaraVicca_old")

CN_SaraVicca3$lon <- round(CN_SaraVicca3$lon,2)
CN_SaraVicca3$lat <- round(CN_SaraVicca3$lat,2)

CN_SaraVicca4 <- aggregate(CN_SaraVicca3, by=list(CN_SaraVicca3$lon,CN_SaraVicca3$lat), mean,na.rm=TRUE)
CN_SaraVicca4 <- CN_SaraVicca4[,c(3,4,5)]

test_new <- merge(NPP_old,CN_SaraVicca4,by=c("lon","lat"),all.x=TRUE) # merging with "old" Sara Vicca's dataset, as far lon and lat both agrees within 0.01 degree.
nrow(test_new)-nrow(NPP_old)

test2 <- subset(test_new,CN_SaraVicca_old>0)
test2 <- test2[,c("site","CN_SaraVicca_old")]
dim(test2)
test2 <- aggregate(test2, by=list(test2$site), mean,na.rm=TRUE)
dim(test2)
test2 <- test2[,c(1,3)]
names(test2) <- c("site","CN_SaraVicca_old")

#now, merge with NPP_final2 - 4 objects need to be merged, at this stage. 
NPP_final2_site <- NPP_final2[,c("lon","lat","z","site","CN_leaf","CN_root")]

NPP_final2_site$no <- 1:nrow(NPP_final2_site)

#merge based on site, for CN_Schulze, CN_Malhi and test2 (old Sara's dataset - for relaxed merged condition - 0.01 deg)
npp_cn1 <-Reduce(function(x,y) merge(x = x, y = y, by = c("site"),all.x=TRUE),
                 list(NPP_final2_site,CN_Schulz2,CN_Malhi2,test2))
nrow(npp_cn1) - nrow(NPP_final2_site)

#merge based on test1 (New Sara's dataset - for strictly merged condition)
npp_cn2 <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat","z"),all.x=TRUE),
                 list(npp_cn1,test1))

npp_cn2 <- npp_cn2[order(npp_cn2$no), ]

#one column have a issue when both test and test2 have data, remove test2
subset(npp_cn2,CN_leaf_sara>0 & CN_SaraVicca_old>0)$no
npp_cn2$CN_SaraVicca_old[npp_cn2$no==subset(npp_cn2,CN_leaf_sara>0 & CN_SaraVicca_old>0)$no] <- NA

summary(npp_cn2$lon - NPP_final2$lon)
summary(npp_cn2)
npp_cn2$CN_stem_final <- npp_cn2$CN_stem_Schulz
npp_cn2$CN_wood_final <- npp_cn2$CN_wood_Schulz

#now, it is the time to combine
npp_cn3 <- npp_cn2 %>% mutate(CN_leaf_final = coalesce(CN_leaf,CN_leaf_Schulz,CN_leaf_alt_malhi,CN_SaraVicca_old,CN_leaf_sara)) %>% # including additional data by secondry C/N merging by Sara
  mutate(CN_leaf_org = coalesce(CN_leaf,CN_leaf_Schulz,CN_leaf_alt_malhi,CN_leaf_sara)) %>% #NOT including additional data by secondry C/N merging by Sara
  mutate(CN_root_final = coalesce(CN_root,CN_root_Schulz))

summary(npp_cn3$CN_leaf_org)
summary(npp_cn3$CN_leaf_final) # it would help us generate 36 more sites

npp_cn4 <- npp_cn3[,c("CN_stem_final","CN_wood_final","CN_leaf_final","CN_leaf_org","CN_root_final")]
hist(npp_cn4$CN_leaf_final)
hist(npp_cn4$CN_stem_final)
hist(npp_cn4$CN_wood_final)
hist(npp_cn4$CN_root_final)

NPP_final3 <- cbind(NPP_final2,npp_cn4)

NPP_final3$lnf_obs_final <-NPP_final3$NPP.foliage/NPP_final3$CN_leaf_final
NPP_final3$lnf_obs_org <-NPP_final3$NPP.foliage/NPP_final3$CN_leaf_org
NPP_final3$bnf_obs_final  <- NPP_final3$BNPP_1/NPP_final3$CN_root_final
NPP_final3$wnf_obs_final  <- NPP_final3$NPP.wood/NPP_final3$CN_wood_final

summary(NPP_final3) 
NPP_final3 <- NPP_final3[,!(names(NPP_final3) %in% c("CN_leaf","CN_root"))]



#add rep_info
#rep_info <- read.csv("~/data/NPP_Yunke/NPP_final_rep.csv")
#NPP_final3$rep_info <- rep_info$rep_info

## 9. add Tiandi's latest forest data (the last data we need!)
tiandi_forest <- read.csv(file="~/data/npp_stoichiometry_forests_tiandi/Site_level_forest_CN_NPP_China_TD_20210104_for_Beni_Yunke.csv")
tiandi_forest <- tiandi_forest[,c(1,2,3,4,11,14,17,18,19,20,21,22)]
head(tiandi_forest)
names(tiandi_forest) <- c("site","lon","lat","z","CN_leaf_final","CN_stem_final","CN_root_final","TNPP_1","ANPP_2",
                          "NPP.foliage","NPP.wood","BNPP_1")

#extend CN_leaf version, as classified in NPP_final3
tiandi_forest$CN_leaf_org <- tiandi_forest$CN_leaf_final

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
tiandi_forest$lnf_obs_org <-tiandi_forest$NPP.foliage/tiandi_forest$CN_leaf_final  
tiandi_forest$bnf_obs_final  <- tiandi_forest$BNPP_1/tiandi_forest$CN_root_final
tiandi_forest$wnf_obs_final  <- tiandi_forest$NPP.wood/tiandi_forest$CN_stem_final #assume stem ratio as wood ratio here?

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

### collect unchanged climate forcing data, based on same lon+lat+z+Begin_year+End_year, 
NPP_Forest <- subset(NPP,pft=="Forest")
dim(NPP_Forest)

NPP_Forest <- NPP_Forest[order(NPP_Forest$site,NPP_Forest$lon,NPP_Forest$lat,NPP_Forest$ANPP_2,NPP_Forest$NPP.foliage,NPP_Forest$Begin_year,NPP_Forest$End_year), ]

#read old one 
NPP_old <- read.csv("/Users/yunpeng/data/forest_npp/NPP_Forest_corrected_Malhi_coord.csv")
NPP_old <- subset(NPP_old,file!="Sara Vicca_stand level")
#NPP_old <- NPP_old[,c("lon","lat","z","Begin_year","End_year","sitename","sitename_fpar")]
NPP_old <- NPP_old[order(NPP_old$site,NPP_old$lon,NPP_old$lat,NPP_old$ANPP_2,NPP_old$NPP.foliage,NPP_old$Begin_year,NPP_old$End_year), ]

summary(NPP_Forest)

summary(NPP_Forest$lon - NPP_old$lon)
summary(NPP_Forest$lat - NPP_old$lat)

#In NPP_foerst Only line 392 and 393 were flipped (see their start_year) - when comparing with NPP_old - let's accept it anyways, since it has no info for rep_info
NPP_Forest[392:393,]
NPP_old[392:393,]
NPP_old$rep_info[392:393]

#in this way, sitename_fpar can be directly merged.
NPP_Forest$sitename_fpar <- NPP_old$sitename_fpar
NPP_Forest$sitename <- NPP_old$sitename
NPP_Forest$rep_info <- NPP_old$rep_info

#ignore below######
#and here is info we need to re-interpolate to update (as recently re-interpolate elevation or measurement year, and using original year info)
corrected_info <- NPP_Forest$sitename[NPP_Forest$z!=NPP_old$z |
                      NPP_Forest$Begin_year!=NPP_old$Begin_year|
                      NPP_Forest$End_year!=NPP_old$End_year]
corrected_info
length(corrected_info) #90
#see this update in "/Users/yunpeng/data/NPP_final/reprocessing_climates/"

#mark it here
NPP_Forest$sitename_replace[NPP_Forest$z!=NPP_old$z |
                              NPP_Forest$Begin_year!=NPP_old$Begin_year|
                              NPP_Forest$End_year!=NPP_old$End_year] <- "yes"

NPP_Forest$file[NPP_Forest$z!=NPP_old$z |
                  NPP_Forest$Begin_year!=NPP_old$Begin_year|
                  NPP_Forest$End_year!=NPP_old$End_year]

#csvfile <- paste("/Users/yunpeng/data/NPP_final/NPP_Forest.csv")
#write_csv(NPP_Forest, path = csvfile)

### Grassland: collect unchanged climate forcing data, based on same lon+lat+z+Begin_year+End_year, 
#old
NPP_grassland  <- read.csv("/Users/yunpeng/data/grassland_npp/NPP_grassland.csv")
#create an old_no to make sure it can be best ordered
NPP_grassland$old_no <- c(1:nrow(NPP_grassland))
NPP_Grassland_old <- NPP_grassland[order(NPP_grassland$lon,NPP_grassland$lat,
                                             NPP_grassland$z,NPP_grassland$ANPP_2,NPP_grassland$TNPP_1), ]

#new
NPP_Grassland_new <- subset(NPP,pft!="Forest")

NPP_Grassland_new <- NPP_Grassland_new[order(NPP_Grassland_new$lon,NPP_Grassland_new$lat,
                                             NPP_Grassland_new$z,NPP_Grassland_new$ANPP_2,NPP_Grassland_new$TNPP_1), ]

summary(NPP_Grassland_new$lon - NPP_Grassland_old$lon)
summary(NPP_Grassland_new$lat - NPP_Grassland_old$lat)
summary(NPP_Grassland_new$z - NPP_Grassland_old$z)
summary(NPP_Grassland_new$ANPP_2 - NPP_Grassland_old$ANPP_2)
summary(NPP_Grassland_new$TNPP_1 - NPP_Grassland_old$TNPP_1)

NPP_Grassland_new$site[NPP_Grassland_new$Begin_year!=NPP_Grassland_old$Begin_year|
                         NPP_Grassland_new$End_year!=NPP_Grassland_old$End_year]
NPP_Grassland_new$Begin_year[NPP_Grassland_new$Begin_year!=NPP_Grassland_old$Begin_year|
                         NPP_Grassland_new$End_year!=NPP_Grassland_old$End_year]
NPP_Grassland_new$End_year[NPP_Grassland_new$Begin_year!=NPP_Grassland_old$Begin_year|
                               NPP_Grassland_new$End_year!=NPP_Grassland_old$End_year]

NPP_Grassland_old$Begin_year[NPP_Grassland_new$Begin_year!=NPP_Grassland_old$Begin_year|
                               NPP_Grassland_new$End_year!=NPP_Grassland_old$End_year]
NPP_Grassland_old$End_year[NPP_Grassland_new$Begin_year!=NPP_Grassland_old$Begin_year|
                             NPP_Grassland_new$End_year!=NPP_Grassland_old$End_year]

#convert them to 2002
NPP_Grassland_new$Begin_year[NPP_Grassland_new$Begin_year!=NPP_Grassland_old$Begin_year|
                               NPP_Grassland_new$End_year!=NPP_Grassland_old$End_year] <- 2002
NPP_Grassland_new$End_year[NPP_Grassland_new$Begin_year!=NPP_Grassland_old$Begin_year|
                             NPP_Grassland_new$End_year!=NPP_Grassland_old$End_year] <- 2002


summary(NPP_Grassland_new$Begin_year - NPP_Grassland_old$Begin_year)
summary(NPP_Grassland_new$End_year - NPP_Grassland_old$End_year)

#in keith's one site- it was set to 2002, while now setting back to 1991-2010, since it is planation - let's ignore it anyways.
#now we can cbind old_no, as originally ranked in old df - it will be used to create sitename later in site simulaion
NPP_Grassland_new$old_no <- NPP_Grassland_old$old_no
NPP_Grassland_new <- NPP_Grassland_new[order(NPP_Grassland_new$old_no), ]
NPP_Grassland_old <- NPP_Grassland_old[order(NPP_Grassland_old$old_no), ]
summary(NPP_Grassland_new$lon - NPP_Grassland_old$lon)
summary(NPP_Grassland_new$lat - NPP_Grassland_old$lat)
summary(NPP_Grassland_new$z - NPP_Grassland_old$z)
summary(NPP_Grassland_new$Begin_year - NPP_Grassland_old$Begin_year)
summary(NPP_Grassland_new$End_year - NPP_Grassland_old$End_year)
summary(NPP_Grassland_new$ANPP_2 - NPP_Grassland_old$ANPP_2)
summary(NPP_Grassland_new$TNPP_1 - NPP_Grassland_old$TNPP_1)

NPP_Grassland_new$rep_info <- NPP_Grassland_old$rep_info
subset(NPP_Grassland_new,rep_info=="rep" | rep_info=="rep2"| rep_info=="rep3")

###########
#now, it is the time to merge with c3c4 information from three different sources
###########
NPP_grassland <- NPP_Grassland_new[,!(names(NPP_Grassland_new) %in% c("age","LAI","observedfAPAR","alpha","soilCN","Evergreen.Deciduous"))]
head(NPP_grassland)

#1. Tian Di's data
tiandi_df_sp <- read.csv("/Users/yunpeng/data/npp_stoichiometry_grasslands_tiandi/China_grassland_CN_stoichiometry_with_matched_NPP_species_legume_20201214.csv")

list_df <- vector(mode = "list", length = nrow(tiandi_df_sp))

for (i in (1:nrow(tiandi_df_sp))){
  list_df[[i]] <- strsplit(tiandi_df_sp$Species_CN[i], "_", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  
}

for (a in (1:nrow(tiandi_df_sp))){
  tiandi_df_sp[a,21:33] <- list_df[[a]][[1]][1:13]
}

t1 <- tiandi_df_sp[,21:33] 
for (i in (1:nrow(t1))){
  t1$no[i] <- i
}

library(reshape)
t2 <- melt(t1, id.vars=c('no'),var='species')
t3 <- na.omit(t2)
t4 <- t3[order(t3$no), ]
t5 <- t4[,c("no","value")]
dim(t5) # number of individuals overall

final_species <- aggregate(no~value,FUN=mean,na.rm=TRUE,data=t5)
dim(final_species) #number of species type

#separate into genus species
for (i in (1:nrow(final_species))){
  final_species[i,3] <- strsplit(final_species$value[i], " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1] #genus
  final_species[i,4] <- strsplit(final_species$value[i], " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2] #species
}
head(final_species)
final_species <- final_species[,c(1,3,4)] 
names(final_species) <- c("speciesname","genus","species")
dim(final_species)
#csvfile <- paste("/Users/yunpeng/data/npp_stoichiometry_grasslands_tiandi/species_name.csv")
#write.csv(final_species, csvfile, row.names = TRUE)

#now, input c3/c4 information from TRY database
c3c4 <- read.csv("/Users/yunpeng/data/c3c4_species/Try20201218143430TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease.csv")
data1 <- c3c4[,c(2,4,5,18)]
dim(data1)
names(data1) <- c("speciesname","genus","species","c3")


final_species2 <- merge(final_species,data1,by=c("speciesname"),all.x=TRUE)

#after having a look at original TRY data, for NA data of final_species2: if the same Genus in TRY database all have recorded c3, then we transfer our NA of same Genus to c3;
# if c3/c4 existed in the same Genus, or Genus is missing, then we set to unknown.
final_species2$c3_final <- final_species2$c3
final_species2$c3_final[6] <- "C4"
final_species2$c3_final[c(16,17,69,70,91,98,100,106,108,110,111,113,114,
                          115,116,117,122,123,124,133,142,153,154,177,178,179,180,206)] <- "unknown"
final_species2$c3_final[final_species2$c3_final==""] <- "tranfered_c3"
final_species2$c3_final[is.na(final_species2$c3_final)==TRUE] <- "tranfered_c3"

#have a look at finalspecies c3c4 data and create a new variable name for percentage
final_species2$c3_percentage <- NA

final_species2$c3_percentage[final_species2$c3_final=="C3"] <- 1
final_species2$c3_percentage[final_species2$c3_final=="tranfered_c3"] <- 1 #with same genus in TRY that all = c3, then also converted to C3
final_species2$c3_percentage[final_species2$c3_final=="C4"] <- 0
final_species2$c3_percentage[final_species2$c3_final=="unknown"] <- NA

#now, time to merge with all individuals data
final_sp_tian <- final_species2[,c("speciesname","c3_percentage")]
names(t5) <- c("no","speciesname")

all_individuals_tian <- merge(t5,final_sp_tian,by=c("speciesname"),all.x=TRUE)
all_individuals_tian <- all_individuals_tian[order(all_individuals_tian$no), ]
#if na occured in certain speciees, then omit (na.rm=TRUE)
all_individuals_tian2 <- aggregate(all_individuals_tian,by=list(all_individuals_tian$no), FUN=mean, na.rm=TRUE)
summary(all_individuals_tian2)
all_individuals_tian2 <- all_individuals_tian2[order(all_individuals_tian2$no), ]

hist(all_individuals_tian2$c3_percentage)

#finally, cbind to original data
tiandi_df_sp2 <- tiandi_df_sp[,1:20]
tiandi_df_sp2$c3_percentage <- all_individuals_tian2$c3_percentage
summary(tiandi_df_sp2)
subset(tiandi_df_sp2,is.na(c3_percentage)==TRUE)

tiandi_df_sp3 <- tiandi_df_sp2[,c("Longitude_CN","Latitude_CN","Altitude_CN","CN_ratio_leaf","ANPP","c3_percentage")]
names(tiandi_df_sp3) <- c("lon","lat","z","CN_leaf_final","ANPP_2","c3_percentage_tiandi")

head(NPP_grassland)

NPP_grassland_final4 <- merge(NPP_grassland,tiandi_df_sp3,by=c("lon","lat","z","CN_leaf_final","ANPP_2"),all.x=TRUE)
summary(NPP_grassland_final4)


#2. Keith's data
keith_c3c4 <- read.csv("/Users/yunpeng/data/NPP_Yunke/NPP_Keith/orig/ABPE.csv")
keith2 <- keith_c3c4[,c("Site","ANPP","C_cycle")]
keith2$c3_percentage_keith[keith2$C_cycle=="C3"] <- 1
keith2$c3_percentage_keith[keith2$C_cycle=="C4"] <- 0
keith2$c3_percentage_keith[keith2$C_cycle=="NA"] <- NA

keith2 <- keith2[,c("Site","ANPP","c3_percentage_keith")]
names(keith2) <- c("site","ANPP_2","c3_percentage_keith")
#merged with site and ANPP_2 (so that each individual is unique).
NPP_grassland_final5 <- merge(NPP_grassland_final4,keith2,by=c("site","ANPP_2"),all.x=TRUE)
summary(NPP_grassland_final5) # 1097-1053 = 44 points were filled now

#3. Campioli
Campioli_c3c4 <- read.csv("/Users/yunpeng/data/campioli/structured_Database1Grassland.csv")
Campioli2 <- Campioli_c3c4[,c("ID","type")]
Campioli2$c3_percentage_Campioli[Campioli2$type=="c3"] <- 1
Campioli2$c3_percentage_Campioli[Campioli2$type=="c4"] <- 0
Campioli2$c3_percentage_Campioli[Campioli2$type=="c3c4"] <- NA # we don't know how much percentage they have (i.e. species number) so we can only extract them from c3c4 percentage map.
Campioli2$c3_percentage_Campioli[Campioli2$type=="NA"] <- NA

Campioli2 <- Campioli2[,c("ID","c3_percentage_Campioli")]
names(Campioli2) <- c("site","c3_percentage_Campioli")
NPP_grassland_final6 <- merge(NPP_grassland_final5,Campioli2,by=c("site"),all.x=TRUE)
summary(NPP_grassland_final6) 
#1097-1012 = 86 points were filled now, which is much less than 142 --> have a look at na table and see which site' missing c3c4 can be filled manually (because the site information given by Cambiopli in two times email are NOT completely the same!!!)
atest <- subset(NPP_grassland_final6,file=="MCampioli" & is.na(c3_percentage_Campioli)==TRUE)
#below is the manually step to fill the c3c4 based on orig c3c4 data given (it was missing in merge because the sitename was not perfectly matached)
NPP_grassland_final6$c3_percentage_Campioli[c(18,19,20,21,22,23,27,28,29,30,31,32,858,859)] <- 1
NPP_grassland_final6$c3_percentage_Campioli[c(807,808,809,810,846,853,856,857,877,878)] <- 0 #US-ccc-D01, US-ccc-D02,US-Kon-D05, US-paw-D01,US-Seg-D01,VE-ori-D01,VE-ori-D02
summary(NPP_grassland_final6) 

##Finally, combine the three above to one dataset
NPP_grassland_final7 <- NPP_grassland_final6 %>% 
  mutate(c3_percentage = coalesce(c3_percentage_tiandi,c3_percentage_keith,c3_percentage_Campioli))
summary(NPP_grassland_final7)

#Alternatively (last), now let's use c3c4 map to fill in those 82 points.
library(rbeni)
c4_still <- as.data.frame(nc_to_df(read_nc_onefile(
  "/Users/yunpeng/data/c4_still/final/c4_percentage.nc"),
  varnam = "c4"))

#use extract function to extract c4 (not gwr!)
names(c4_still) <- c("lon","lat","c4")
coordinates(c4_still) <- ~lon+lat 
gridded(c4_still) <- TRUE
rc4_global <- raster(c4_still, "c4") 

#aggregate based on lon and lat firstly
NPP_grassland_final7_site <- aggregate(NPP_grassland_final7,by=list(NPP_grassland_final7$lon,NPP_grassland_final7$lat), FUN=mean, na.rm=TRUE) #site-mean
NPP_grassland_final7_site <- NPP_grassland_final7_site[,c("lon","lat")]

sp_sites <- SpatialPoints(NPP_grassland_final7_site) # only select lon and lat

NPP_grassland_final7_site <- raster::extract(rc4_global, sp_sites, sp = TRUE) %>% as_tibble() %>% 
  right_join(NPP_grassland_final7_site, by = c("lon", "lat")) %>% 
  dplyr::rename( c4_percentage_map = c4)
dim(NPP_grassland_final7_site)
hist(NPP_grassland_final7_site$c4_percentage_map) # most c4 percentage =0, which is great.

#now, merge back to site
NPP_grassland_final7_site$c3_precentage_map <- 1-NPP_grassland_final7_site$c4_percentage_map
NPP_grassland_final7_site <- NPP_grassland_final7_site[,c("lon","lat","c3_precentage_map")]

NPP_grassland_final8 <- merge(NPP_grassland_final7,NPP_grassland_final7_site,by=c("lon","lat"),all.x=TRUE)
dim(NPP_grassland_final8)

summary(NPP_grassland_final8)

#compare measured and predicted c3 percentage --> very different!
plot(NPP_grassland_final8$c3_percentage~NPP_grassland_final8$c3_precentage_map)

NPP_grassland_final8 <- NPP_grassland_final8[order(NPP_grassland_final8$old_no), ]

#now, combine them: primary based on measured data, alternatively based on map data.
NPP_grassland_final9 <- NPP_grassland_final8 %>% 
  mutate(c3_percentage_final = coalesce(c3_percentage,c3_precentage_map))
#which points were filled by map
subset(NPP_grassland_final9,is.na(c3_percentage)==TRUE)$c3_percentage_final

#additional step: 
#add more gpp data from Compioli et al. SI table 1, and remove repeated data
#NPP_grassland_final9_v2
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="IT-bea-D02"] <- 1568
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-che-D01"] <- 626
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="DE-gri-D01"] <- 1233#repeated
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="CN-Hab-F01"] <- 634
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="RU-ha1-F01"] <- 519#repeated
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="RU-ha3-F01"] <- 526 #repeated
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="CN-Inn-D01_C"] <- 182
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-kbs-D01"] <- 1015
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-kbs-D04"] <- 512
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-kbs-D05"] <- 374
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-kbs-D03"] <- 793
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-jas-D01"] <- 516
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-kon-D05"] <- 1151
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="RU-krs-D01"] <- 1611
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="CA-Let-F01"] <- 280
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="CA-mat-D01"] <- 786
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-osg-D01"] <- 1890
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="CG-tch-D01"] <- 1572
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="US-Spe-D01"] <- 829

#interpolate gpp to Cambioli's data from keith's source (with same site name)
NPP_grassland_final9$GPP[NPP_grassland_final9$site=="CN-Inn-F01"] <- 204 #this value is reasonable comparing with NPP/GPP
#NPP_grassland_final9$GPP[NPP_grassland_final9$site=="KZ-shr-D01"] <- 357 #incosistent between keith and cambiopli's data! (i.e. leading to unreasonable values with GPP < TNPP)
#NPP_grassland_final9$GPP[NPP_grassland_final9$site=="RU-ha2-F01"] <- 648 #incosistent between keith and cambiopli's data!  (i.e. leading to unreasonable values with GPP slightly higher than TNPP - CUE =94%)

NPP_grassland_final9$BNPP_1 <- NPP_grassland_final9$TNPP_1 - NPP_grassland_final9$ANPP_2
NPP_grassland_final10 <- NPP_grassland_final9

#add more data from /Users/yunpeng/data/NPP_Yunke/NPP_Keith/orig/Bremer and Ham 2010.pdf from keith's grassland data
#based on its Table 3 - GPP, for example, at the site below, was calculated as averages of 2 values (at 2 period) from site BA, so does site BB. 
NPP_grassland_final10$TNPP_1[NPP_grassland_final10$site=="US-Kon-D02"] <- ((1669-1354) + (2269-1666))/2
NPP_grassland_final10$TNPP_1[NPP_grassland_final10$site=="US-Kon-D03"] <- ((1368-1185) + (1997-1495))/2

NPP_grassland_final10$BNPP_1[NPP_grassland_final10$site=="US-Kon-D02"] <- NPP_grassland_final10$TNPP_1[NPP_grassland_final10$site=="US-Kon-D02"] - NPP_grassland_final10$ANPP_2[NPP_grassland_final10$site=="US-Kon-D02"]
NPP_grassland_final10$BNPP_1[NPP_grassland_final10$site=="US-Kon-D03"] <- NPP_grassland_final10$TNPP_1[NPP_grassland_final10$site=="US-Kon-D03"] - NPP_grassland_final10$ANPP_2[NPP_grassland_final10$site=="US-Kon-D03"]

summary(NPP_grassland_final10)

NPP_grassland_final11 <- NPP_grassland_final10
NPP_grassland_final12 <- NPP_grassland_final11[,!(names(NPP_grassland_final11) %in% c("c3_percentage_tiandi","c3_percentage_keith","c3_percentage_Campioli","c3_precentage_map"))]
NPP_grassland_final12 <- NPP_grassland_final12[order(NPP_grassland_final12$old_no), ]

#double check if two df is the same - to cbind with sitename later on in site simulation in another file: grassland_simulation
NPP_grassland  <- read.csv("/Users/yunpeng/data/grassland_npp/NPP_grassland.csv")
summary(NPP_grassland_final12$lon - NPP_grassland$lon)
summary(NPP_grassland_final12$lat - NPP_grassland$lat)
summary(NPP_grassland_final12$z - NPP_grassland$z)
summary(NPP_grassland_final12$Begin_year - NPP_grassland$Begin_year)
summary(NPP_grassland_final12$End_year - NPP_grassland$End_year)
summary(NPP_grassland_final12$ANPP_2 - NPP_grassland$ANPP_2)
summary(NPP_grassland_final12$TNPP_1 - NPP_grassland$TNPP_1)

#csvfile <- paste("/Users/yunpeng/data/NPP_Grassland_final/NPP_grassland.csv")
#write_csv(NPP_grassland_final12, path = csvfile)


#now, newly add net minerlization rate and GCME N uptake sites
#### Input N uptake
#(1) newly added Nmin rate data from Finzi
Finzi <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nmin_Finzi/Nmin_Finzi.csv")
names(Finzi)[names(Finzi) == "Lat"] <- "lat"
names(Finzi)[names(Finzi) == "Long"] <- "lon"
devtools::load_all("/Users/yunpeng/yunkepeng/Grassland_new_ingestr_rsofun_20210326/ingestr/")

Finzi$pft <- "Forest"
Finzi$pft[Finzi$Biome=="temp grass"]<- "Grassland"

Finzi_Forest_sitemean <- aggregate(Finzi,by=list(Finzi$lon,Finzi$lat), FUN=mean, na.rm=TRUE) #site-mean
dim(Finzi_Forest_sitemean)
for (i in 1:nrow(Finzi_Forest_sitemean)){
  Finzi_Forest_sitemean$sitename[i] <- paste("Finzi_Forest",i,sep = "") # this is also sitename for fpar
}
df_etopo <- ingest(Finzi_Forest_sitemean,source = "etopo1",dir = "~/data/etopo/" )
Finzi_Forest_sitemean$z <- as.numeric(as.data.frame(df_etopo$data))
Finzi_Forest_sitemean$z[Finzi_Forest_sitemean$z< 0] <- 0

Finzi_Forest_sitemean2 <- Finzi_Forest_sitemean[,c("lon","lat","z","sitename")]

Finzi_all <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat"),all.x=TRUE), 
                   list(Finzi,Finzi_Forest_sitemean2))

Finzi_all$Begin_year <- 1984
Finzi_all$End_year <- 2013

Finzi_all <- Finzi_all[,c("lon","lat","z","pft","Nmin","Begin_year","End_year","sitename")]
names(Finzi_all) <- c("lon","lat","z","pft","Nuptake","Begin_year","End_year","site")
Finzi_all$file <- "Finzi, Nmin"

head(Finzi_all)

#(2) Nuptake from gcme
gcme_nuptake <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nuptake_gcme/gcme_nuptake_coord_interpolated.csv")
for (i in 1:nrow(gcme_nuptake)){
  gcme_nuptake$sitename[i] <- paste("gcme",i,sep = "")
}

#elv
df_etopo <- ingest(gcme_nuptake,source = "etopo1",dir = "~/data/etopo/" )
gcme_nuptake$z <- as.numeric(as.data.frame(df_etopo$data))
gcme_nuptake$z[gcme_nuptake$z<0] <- 0

siteinfo_gcme <- data.frame(
  sitename = gcme_nuptake$sitename,
  lon = gcme_nuptake$lon,
  lat = gcme_nuptake$lat,
  elv = gcme_nuptake$z,
  year_start = gcme_nuptake$start_yr,
  year_end = gcme_nuptake$start_yr + gcme_nuptake$Year_long -1
)

siteinfo_gcme$exp_nam <- gcme_nuptake$exp_nam

gcme_data <- read.csv("/Users/yunpeng/data/NPP_Yunke/Nuptake_gcme/gcme_nuptake_data.csv")
gcme_data <- gcme_data[,c("exp_nam","ambient","Unit")]
#all converting to gN/m2/yr
gcme_data$ambient[gcme_data$Unit=="Kg_N_ha-1"] <- gcme_data$ambient[gcme_data$Unit=="Kg_N_ha-1"]/10
gcme_data$ambient[gcme_data$Unit=="kg_N/ha"] <- gcme_data$ambient[gcme_data$Unit=="kg_N/ha"]/10
gcme_data$ambient[gcme_data$Unit=="mg_N/kg*day_"] <- NA #quite weired about the unit, for one site. Disregard them first
hist(gcme_data$ambient)
#why some values are too low?
gcme_data <- subset(gcme_data,ambient>0)
subset(gcme_data,ambient<1) %>% group_by(exp_nam) %>% summarise(number = n())
#after look, "RiceFACE_Japan_A_1998_39,40_141" looks fine, as most of them were still in good range. But remove the other 5 sites (they may be wrong due to measurement mistakes or unit errors, we don't know)
gcme_data_final <- subset(gcme_data,exp_nam!="Michigan_UNDERC_bog" & exp_nam!="Michigan_UNDERC_intermFen" & exp_nam!="Michigan_UNDERC_richFen"&
                            exp_nam!="RiceFACE_China_32N_120E_Or_Tr_7"& exp_nam!="TL_7")
hist(gcme_data_final$ambient)

#finally, mergeing them to obtain siteinfo
gcme_data_final_forest <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam"),all.x=TRUE),list(gcme_data_final,siteinfo_gcme))

gcme_data_final_forest <- gcme_data_final_forest[,c("exp_nam","ambient","lon","lat","elv","year_start","year_end")]
names(gcme_data_final_forest) <- c("site","Nuptake","lon","lat","z","Begin_year","End_year")
gcme_data_final_forest$file <- "GCME, Nuptake"
gcme_data_final_forest$pft <- "Forest" 
#rbind them to get final Nmin. data
Nmin_final <- dplyr::bind_rows(Finzi_all,gcme_data_final_forest)
Nmin_final

#Final dataset to be ingested together
dim(NPP_Forest)
dim(NPP_grassland_final12)
dim(Nmin_final)

siteinfo_final <- dplyr::bind_rows(NPP_Forest,NPP_grassland_final12,Nmin_final)

siteinfo_final_ingest <- aggregate(siteinfo_final,by=list(siteinfo_final$lon,siteinfo_final$lat,siteinfo_final$z,siteinfo_final$Begin_year,siteinfo_final$End_year), FUN=mean, na.rm=TRUE)
siteinfo_final_ingest <- siteinfo_final_ingest[,c("lon","lat","z","Begin_year","End_year")]

siteinfo_final_ingest$year_start <- siteinfo_final_ingest$Begin_year
siteinfo_final_ingest$year_end <- siteinfo_final_ingest$End_year
siteinfo_final_ingest$year_start[siteinfo_final_ingest$Begin_year<1979] <- 1979
siteinfo_final_ingest$year_end[siteinfo_final_ingest$Begin_year<1979] <- 1988

summary(siteinfo_final_ingest)

for (i in 1:nrow(siteinfo_final_ingest)){
  siteinfo_final_ingest$sitename_new[i] <- paste("final",i,sep = "") # this is also sitename for fpar
}
dim(siteinfo_final_ingest)

siteinfo_final_ingest <-  siteinfo_final_ingest %>% dplyr::mutate(date_start = lubridate::ymd(paste0(year_start, "-01-01"))) %>%
  dplyr::mutate(date_end = lubridate::ymd(paste0(year_end, "-12-31"))) 

head(siteinfo_final_ingest)

aaa <-Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat","z","Begin_year","End_year"),all.x=TRUE), 
                   list(siteinfo_final,siteinfo_final_ingest[,c("lon","lat","z","Begin_year","End_year","sitename_new")]))

