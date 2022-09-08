library(lme4)
library(MuMIn)
library(lmerTest)
library(visreg)
library(maps)
library(dplyr)
library(dplyr)
library(maps)
library(rworldmap)
library(readr)
library(lme4)
library(MuMIn)
library(lmerTest)
library(visreg)
library(maps)
library(dplyr)
library(cowplot)
library(ggplot2)
#Input climate data ("input site_information.csv")
devtools::load_all("/Users/yunpeng/yunkepeng/latest_packages/rbeni/") # using beni's latest package.

#### Input all indiduals data, and applied pre-processing ("input individuals.csv")
SP_input <- read.csv("/Users/yunpeng/data/leaf_traits/combined_leaf_traits_updated.csv")

#SP_input <- subset(SP_input,source!="Bahar et al 2017 New Phytologist")
#if remove Bahar - set soil coordinates resolution as 1, so that to get more points and not being filtered (atkin and bahar points are not having same coordinates that cannot be consistenly merged with soil - so annoyed!).
#in this way, Ptotal - leaf N is still signifciant in site-species data
#if not removing (like now) - set soil coordinates resolution as 2

SP_input$z[SP_input$z<0] <- 0
SP_input3 <- SP_input

#Read soil data - alraedy tested this with old collection - highly correlated - but this method is more precise. 
#1.read Lloyd
soil_lloyd <- read.csv("/Users/yunpeng/data/soil/Lloyd/Lloyd.csv")
soil_lloyd <- soil_lloyd[,c("Plot","lon","lat","pH","P_total","CN")]

#2. read TROBIT
soil1_trobit <-read.csv("/Users/yunpeng/data/soil/TROBIT/trobit_data_for_ICP.csv")
soil1_trobit <- soil1_trobit[,c("Plot","CN","pH","P_total")]

soil2_trobit <- read.csv("/Users/yunpeng/data/leaf_traits/TROBIT/orig/TROBIT shared data.csv")
soil2_trobit <- soil2_trobit[,c("site","lat","lon")]
names(soil2_trobit)<- c("Plot","lat","lon")
soil3_trobit <-Reduce(function(x,y) merge(x = x, y = y, by = c("Plot"),all.x=TRUE), 
                      list(soil1_trobit,soil2_trobit))
subset(soil3_trobit,is.na(lat)==TRUE)
#correct such lat/lon manually - based on some records from /Users/yunpeng/data/soil/allsoil.csv and /Users/yunpeng/data/soil/Lloyd/Lloyd.csv
soil3_trobit$lat[soil3_trobit$Plot=="BDA-01"] <- 10.939685; soil3_trobit$lon[soil3_trobit$Plot=="BDA-01"] <- -3.149486
soil3_trobit$lat[soil3_trobit$Plot=="BDA-02"] <- 10.939893; soil3_trobit$lon[soil3_trobit$Plot=="BDA-02"] <- -3.15431
soil3_trobit$lat[soil3_trobit$Plot=="BDA-03"] <- 10.86532; soil3_trobit$lon[soil3_trobit$Plot=="BDA-03"] <- -3.072615
soil3_trobit$lat[soil3_trobit$Plot=="SIN-01"] <- NA; soil3_trobit$lon[soil3_trobit$Plot=="SIN-01"] <- NA
soil3_trobit$lat[soil3_trobit$Plot=="TAN-01"] <- NA; soil3_trobit$lon[soil3_trobit$Plot=="TAN-01"] <- NA
soil3_trobit$lat[soil3_trobit$Plot=="TAP-123"] <- -3.309; soil3_trobit$lon[soil3_trobit$Plot=="TAP-123"] <- -54.94

soil_final <- dplyr::bind_rows(soil3_trobit, soil_lloyd)

#aggregate based on lon/lat
soil_final$lon_2<- round(soil_final$lon,2)
soil_final$lat_2<- round(soil_final$lat,2)

soil_final_sitemean <- aggregate(soil_final,by=list(soil_final$lon_2,soil_final$lat_2), FUN=mean, na.rm=TRUE) #site-mean
dim(soil_final_sitemean)
soil_final_sitemean <- soil_final_sitemean[,c("lon_2","lat_2","CN","pH","P_total")]

SP_input3$lon_2<- round(SP_input3$lon,2)
SP_input3$lat_2<- round(SP_input3$lat,2)

SP_input4 <- merge(SP_input3,soil_final_sitemean,by=c("lon_2","lat_2"),all.x=TRUE) #merged sitename to SP data
dim(SP_input3)
dim(SP_input4)
SP_input4$source[SP_input4$source=="Dong Ning collection"] <- "Wright"

summary(SP_input4)

SP_input4$nmass <- SP_input4$narea/SP_input4$lma
SP_input4$leafCN <- SP_input4$C_percent/100/SP_input4$nmass 
SP_input4$cmass <- SP_input4$C_percent/100
SP_input4$carea <- SP_input4$cmass/SP_input4$lma

SP_input4_sitemean <- aggregate(SP_input4,by=list(SP_input4$lon_2,SP_input4$lat_2), FUN=mean, na.rm=TRUE) #site-mean

larger_size <- theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                     plot.subtitle=element_text(size=15))

a1 <- analyse_modobs2(SP_input4_sitemean,"CN","narea", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size
summary(lm(narea~CN,SP_input4_sitemean))

a2 <- analyse_modobs2(SP_input4_sitemean,"CN","nmass", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(nmass~CN,SP_input4_sitemean))

a3 <- analyse_modobs2(SP_input4_sitemean,"CN","leafCN", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(leafCN~CN,SP_input4_sitemean))

a4 <- analyse_modobs2(SP_input4_sitemean,"CN","lma", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(lma~CN,SP_input4_sitemean))

a5 <- analyse_modobs2(SP_input4_sitemean,"CN","Vcmax25", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size
summary(lm(Vcmax25~CN,SP_input4_sitemean))

a6 <- analyse_modobs2(SP_input4_sitemean,"CN","cmass", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(cmass~CN,SP_input4_sitemean))

a7 <- analyse_modobs2(SP_input4_sitemean,"CN","carea", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(carea~CN,SP_input4_sitemean))

b1 <- analyse_modobs2(SP_input4_sitemean,"P_total","narea", type = "points",relative=TRUE)$gg+labs(x ="Soil P")+larger_size
summary(lm(narea~P_total,SP_input4_sitemean))

b2 <- analyse_modobs2(SP_input4_sitemean,"P_total","nmass", type = "points",relative=TRUE)$gg+labs(x ="Soil P")+larger_size
summary(lm(nmass~P_total,SP_input4_sitemean))

b3 <- analyse_modobs2(SP_input4_sitemean,"P_total","leafCN", type = "points",relative=TRUE)$gg+labs(x ="Soil P")+larger_size
summary(lm(leafCN~P_total,SP_input4_sitemean))

b4 <- analyse_modobs2(SP_input4_sitemean,"P_total","lma", type = "points",relative=TRUE)$gg+labs(x ="Soil P")+larger_size
summary(lm(lma~P_total,SP_input4_sitemean))

b5 <- analyse_modobs2(SP_input4_sitemean,"P_total","Vcmax25", type = "points",relative=TRUE)$gg+labs(x ="Soil P")+larger_size
summary(lm(Vcmax25~P_total,SP_input4_sitemean))

b6 <- analyse_modobs2(SP_input4_sitemean,"P_total","cmass", type = "points",relative=TRUE)$gg+labs(x ="Soil P")+larger_size
summary(lm(cmass~P_total,SP_input4_sitemean))

b7 <- analyse_modobs2(SP_input4_sitemean,"P_total","carea", type = "points",relative=TRUE)$gg+labs(x ="Soil P")+larger_size
summary(lm(carea~P_total,SP_input4_sitemean))

c1 <- analyse_modobs2(SP_input4_sitemean,"pH","narea", type = "points",relative=TRUE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(narea~pH,SP_input4_sitemean))

c2 <- analyse_modobs2(SP_input4_sitemean,"pH","nmass", type = "points",relative=TRUE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(nmass~pH,SP_input4_sitemean))

c3 <- analyse_modobs2(SP_input4_sitemean,"pH","leafCN", type = "points",relative=TRUE)$gg+labs(x ="pH")+larger_size
summary(lm(leafCN~pH,SP_input4_sitemean))

c4 <- analyse_modobs2(SP_input4_sitemean,"pH","lma", type = "points",relative=TRUE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(lma~pH,SP_input4_sitemean))

c5 <- analyse_modobs2(SP_input4_sitemean,"pH","Vcmax25", type = "points",relative=TRUE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(Vcmax25~pH,SP_input4_sitemean))

c6 <- analyse_modobs2(SP_input4_sitemean,"pH","cmass", type = "points",relative=TRUE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(cmass~pH,SP_input4_sitemean))

c7 <- analyse_modobs2(SP_input4_sitemean,"pH","carea", type = "points",relative=TRUE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)
summary(lm(carea~pH,SP_input4_sitemean))

plot_grid(a1,a2,a3,a4,a5,a6,a7,
          b1,b2,b3,b4,b5,b6,b7,
          c1,c2,c3,c4,c5,c6,c7,
          nrow=3)
ggsave(paste("~/data/soilN_P.jpg",sep=""),width = 30, height = 15)

summary(lm(CN~pH,SP_input4_sitemean))

#check site info
unique(subset(SP_input4,CN>0 & leafCN>0)$site)
unique(subset(SP_input4,CN>0 & leafCN>0)$source)

unique(subset(SP_input4,CN>0 & narea>0)$site)
unique(subset(SP_input4,CN>0 & narea>0)$source)
#bahar, atkin, bloomfield and TROBIT may lack leaf C data?


#check 1
summary(lm(leafCN~log(CN),SP_input4_sitemean))

SP_input4_sitemean$vcmax25_lma <- SP_input4_sitemean$Vcmax25/SP_input4_sitemean$lma
summary(lm(nmass~CN+vcmax25_lma,SP_input4_sitemean))
