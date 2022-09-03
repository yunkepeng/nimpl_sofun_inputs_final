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
SP_input4_sitemean <- aggregate(SP_input4,by=list(SP_input4$lon_2,SP_input4$lat_2), FUN=mean, na.rm=TRUE) #site-mean


larger_size <- theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                     plot.subtitle=element_text(size=15))

a1 <- analyse_modobs2(SP_input4_sitemean,"CN","narea", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size
summary(lm(narea~CN,SP_input4_sitemean))

a2 <- analyse_modobs2(SP_input4_sitemean,"CN","nmass", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size
summary(lm(nmass~CN,SP_input4_sitemean))

a3 <- analyse_modobs2(SP_input4_sitemean,"CN","leafCN", type = "points",relative=TRUE)$gg+labs(x ="Soil C/N")+larger_size
summary(lm(leafCN~CN,SP_input4_sitemean))

plot_grid(a1,a2,a3, nrow=1)
ggsave(paste("~/data/soilCN.jpg",sep=""),width = 15, height = 5)

#check site info
unique(subset(SP_input4,CN>0 & leafCN>0)$site)
unique(subset(SP_input4,CN>0 & leafCN>0)$source)

unique(subset(SP_input4,CN>0 & narea>0)$site)
unique(subset(SP_input4,CN>0 & narea>0)$source)
#bahar, atkin, bloomfield and TROBIT may lack leaf C data?