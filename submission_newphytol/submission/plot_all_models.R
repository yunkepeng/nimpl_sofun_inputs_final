# Setup 
library(MLmetrics)
library(vhs) ##remotes::install_github("cj-holmes/vhs")
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
library(caret)
library(recipes)
devtools::load_all("/Users/yunpeng/yunkepeng/latest_packages/rbeni/") # using beni's latest package.
#library(rbeni)
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)
library(lme4)
library(lmerTest)
library("PerformanceAnalytics")
library(MuMIn)
library(tidyverse)
library(ggplot2)
library(lme4)
library(visreg)
library(ggpubr)
library(car)
library("ggplotify")
library(remotes)
library(tune)
library(relaimpo)

source("/Users/yunpeng/yunkepeng/CNuptake_MS/R/cal_nue.R")
source("/Users/yunpeng/yunkepeng/CNuptake_MS/R/analyse_modobs2.R")
source("/Users/yunpeng/yunkepeng/CNuptake_MS/R/stepwise.R")
source("/Users/yunpeng/yunkepeng/CNuptake_MS/R/stepwise_lm.R")
source("/Users/yunpeng/yunkepeng/CNuptake_MS/R/calc_area.R")

# reset validation metrics info
white <- theme(plot.background=element_rect(fill="white", color="white"))

larger_size <- theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                     plot.subtitle=element_text(size=15))

# Empirical models

## Read data and pre-processing


#remove alpha throughout all analyses because 
#(1) D and alpha shows contridictory results for NUE 
#(2) D and alpha repeated in model selection

#input NPP and nmin file and combine
npp_dataset <- read.csv("~/data/NPP_Yunke/NPP_dataset.csv")
nmin_dataset <- read.csv("~/data/Nmin_Finzi/Nmin_dataset.csv")
NPP_all <- dplyr::bind_rows(npp_dataset, nmin_dataset) 

#summarise number of sites
dim(subset(NPP_all,is.na(Nmin)==TRUE) %>% group_by(site)  %>% summarise(mean = mean(lon)))

NPP_all$NPP.foliage[NPP_all$NPP.foliage==0] <-NA
NPP_all$NPP.wood[NPP_all$NPP.wood==0] <-NA
NPP_all$ANPP_2[NPP_all$ANPP_2==0] <-NA
NPP_all$age[NPP_all$age==0] <-NA

NPP_all$site_a <- NPP_all$site
NPP_all$tnpp_a <- NPP_all$TNPP_1
NPP_all$anpp_a <- NPP_all$ANPP_2

NPP_all$anpp_tnpp_a <-  log((NPP_all$ANPP_2/NPP_all$tnpp_a)/(1-(NPP_all$ANPP_2/NPP_all$tnpp_a)))
NPP_all$anpp_leafnpp_a <-  log((NPP_all$NPP.foliage/NPP_all$ANPP_2)/(1-(NPP_all$NPP.foliage/NPP_all$ANPP_2)))

NPP_all$soilCN_a <- log(NPP_all$soilCN)
NPP_all$observedfAPAR_a <- NPP_all$observedfAPAR
NPP_all$obs_age_a <- log(NPP_all$age)

NPP_all$age_a <- log(NPP_all$mapped_age)
NPP_all$Tg_a <- NPP_all$Tg
NPP_all$PPFD_a <- log(NPP_all$PPFD)
NPP_all$vpd_a <- log(NPP_all$vpd)
NPP_all$fAPAR_a <- NPP_all$fAPAR
NPP_all$CNrt_a <- log(NPP_all$CNrt)
NPP_all$LMA_a <- log(NPP_all$LMA)
NPP_all$vcmax25_a <- log(NPP_all$vcmax25)
NPP_all$ndep_a <- log(NPP_all$ndep)

NPP_forest <- subset(NPP_all,pft=="Forest")


#set leaf Cmass, wood C/N, root C/N and wood tissue percentage in forest as constant
#see below code
leaf_c_forest <- 0.47 

#wood_cn_forest <- 100
#root_cn_forest <- 94
#wood_tissue_percentage <- 1

###Below was not used
#set wood and root ratio (from mean values of Zhang et al. 2019 GCB doi: 10.1111/gcb.14973)
wood_cn_forest <- 319.04
root_cn_forest <- 94
#set mean values of living tissue percentage of wood in Morris et al (2016) New Phytologist
# mean value of ray and axial parenchyma (RAP) = 27.2%
wood_tissue_percentage <- 1 #0.272
###

#For 18 and 41 they are fixed values of leaf and root C/N, as derived from "Tian Di Grassland"
#see line 884 and 885 in preparation code: https://github.com/yunkepeng/nimpl_sofun_inputs_final/blob/main/submission_newphytol/preprocessing/Forest_site_orig.R, 
#data all derived from "TianDi Grassland", but their BP, ANPP, BNPP data not used in our study
#summary(aggregate(subset(dataset6,pft=="Grassland"),by=list(subset(dataset6,pft=="Grassland")$site), FUN=mean, na.rm=TRUE)$CN_leaf_final)
#summary(aggregate(subset(dataset6,pft=="Grassland"),by=list(subset(dataset6,pft=="Grassland")$site), FUN=mean, na.rm=TRUE)$CN_root_final)
leaf_cn_grassland <- 18

root_cn_grassland <- 41

#nre constant at grassland - see below code - median values of literature Du and Deng
nre_grassland_constant <- 0.69

#applying null model of anpp/bp - see below code
anpp_bp_grassland <- 0.50


### Map with sites
#check why some grassland BP and ANPP is so high
outliers <- subset(NPP_all,pft=="Grassland" & ANPP_2>500)
newmap <- getMap(resolution = "low")
sp::plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(outliers$lon,outliers$lat, col="green", pch=16,cex=2)



## Fit NPP models

### BP
BP_dataset <- na.omit(NPP_forest[,c("tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")])
#model1 <- stepwise(BP_dataset,"tnpp_a")
#model1[[1]]
#model1[[2]]
#bp_model <- (lmer(tnpp_a~Tg_a+observedfAPAR_a+obs_age_a+PPFD_a+alpha_a+(1|site_a),data=BP_dataset))
#summary(bp_model)

#ndep check - BP adding ndep - significant
BP_dataset_ndep <- na.omit(NPP_forest[,c("tnpp_a","ndep_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(BP_dataset_ndep)
stepwise(BP_dataset_ndep,"tnpp_a")[[1]]
stepwise(BP_dataset_ndep,"tnpp_a")[[3]]
bp_model_ndep <- (lmer(tnpp_a~Tg_a+fAPAR_a+ndep_a+CNrt_a+age_a+(1|site_a),data=BP_dataset_ndep))
summary(bp_model_ndep)
r.squaredGLMM(bp_model_ndep)

#latitude check - check if it can replace Tg as a good predictor?
NPP_forest$lat_abs <- abs(NPP_forest$lat)
BP_dataset_lat <- na.omit(NPP_forest[,c("tnpp_a","lat_abs","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(BP_dataset_lat)
stepwise(BP_dataset_lat,"tnpp_a")[[1]]
stepwise(BP_dataset_lat,"tnpp_a")[[3]]
#lat will not help us to improve the model!


### ANPP
#ndep check - ANPP/BP adding ndep - works!
anpp_tnpp_dataset_ndep <- na.omit(NPP_forest[,c("anpp_tnpp_a","ndep_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
stepwise(anpp_tnpp_dataset_ndep,"anpp_tnpp_a")[[1]]
anpp_tnpp_model_ndep <- (lmer(anpp_tnpp_a~ndep_a+CNrt_a+PPFD_a+Tg_a+(1|site_a),data=anpp_tnpp_dataset_ndep))
summary(anpp_tnpp_model_ndep)
r.squaredGLMM(anpp_tnpp_model_ndep)
#ndep check - leaf.npp/ANPP adding ndep - not improved! Ndep is non-significant
anpp_leafnpp_dataset_noage_ndep <- na.omit(NPP_forest[,c("anpp_leafnpp_a","ndep_a","Tg_a","PPFD_a","vpd_a","site_a")])
stepwise(anpp_leafnpp_dataset_noage_ndep,"anpp_leafnpp_a")[[1]]
r.squaredGLMM(lmer(anpp_leafnpp_a~Tg_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_noage_ndep))
summary(lmer(anpp_leafnpp_a~ndep_a+vpd_a+Tg_a+(1|site_a),data=anpp_leafnpp_dataset_noage_ndep))
r.squaredGLMM(lmer(anpp_leafnpp_a~ndep_a+vpd_a+Tg_a+(1|site_a),data=anpp_leafnpp_dataset_noage_ndep))

### BP vs. model predictors

#now, start works
BP_dataset2 <- na.omit(NPP_forest[,c("tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(BP_dataset2)
a2 <- stepwise(BP_dataset2,"tnpp_a")
a2[[1]]
a2[[2]]
a2[[3]]
bp_model <- (lmer(tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+(1|site_a),data=BP_dataset2))
summary(bp_model)
r.squaredGLMM(bp_model)

#check how many data were removed
nrow(BP_dataset)/nrow(BP_dataset2)

vif_bp <- vif((lmer(tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=BP_dataset2)))

#check significance at iteration
summary((lmer(tnpp_a~Tg_a+(1|site_a),data=BP_dataset2)))
summary((lmer(tnpp_a~Tg_a+fAPAR_a+(1|site_a),data=BP_dataset2)))
summary((lmer(tnpp_a~Tg_a+fAPAR_a+age_a+(1|site_a),data=BP_dataset2)))
summary((lmer(tnpp_a~Tg_a+fAPAR_a+age_a+CNrt_a+(1|site_a),data=BP_dataset2)))
summary((lmer(tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+(1|site_a),data=BP_dataset2)))
summary((lmer(tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=BP_dataset2)))

### ANPP vs. modelled predictors

anpp_tnpp_dataset2 <- na.omit(NPP_forest[,c("anpp_tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(subset(NPP_forest,ANPP_2>0))
model2a <- stepwise(anpp_tnpp_dataset2,"anpp_tnpp_a")
model2a[[1]]
model2a[[2]]
model2a[[3]]
anpp_tnpp_model <- (lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+age_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(anpp_tnpp_model)
r.squaredGLMM(anpp_tnpp_model)

vif_anpp_tnpp <- vif((lmer(anpp_tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=anpp_tnpp_dataset2)))

#check significance at iteration
summary(lmer(anpp_tnpp_a~CNrt_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(lmer(anpp_tnpp_a~CNrt_a+PPFD_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+age_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+age_a+vpd_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+age_a+vpd_a+fAPAR_a+(1|site_a),data=anpp_tnpp_dataset2))


### Leaf NPP
#with age, but age should be removed since it shows higher AIC and lower R2
anpp_leafnpp_dataset_age <- na.omit(NPP_forest[,c("anpp_leafnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
model3a <- stepwise(anpp_leafnpp_dataset_age,"anpp_leafnpp_a")
model3a[[1]]
model3a[[2]]
model3a[[3]]

#check significance at iteration
summary(lmer(anpp_leafnpp_a~age_a+(1|site_a),data=anpp_leafnpp_dataset_age))
summary(lmer(anpp_leafnpp_a~age_a+PPFD_a+(1|site_a),data=anpp_leafnpp_dataset_age))
summary(lmer(anpp_leafnpp_a~age_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 
summary(lmer(anpp_leafnpp_a~age_a+PPFD_a+vpd_a+Tg_a+(1|site_a),data=anpp_leafnpp_dataset_age)) # non-sigificant in age, so age should be removed
summary(lmer(anpp_leafnpp_a~age_a+PPFD_a+vpd_a+Tg_a+fAPAR_a+(1|site_a),data=anpp_leafnpp_dataset_age))
summary(lmer(anpp_leafnpp_a~age_a+PPFD_a+vpd_a+Tg_a+fAPAR_a+CNrt_a+(1|site_a),data=anpp_leafnpp_dataset_age))


#without age - re-selection - this is best
anpp_leafnpp_dataset <- na.omit(NPP_forest[,c("anpp_leafnpp_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
model3 <- stepwise(anpp_leafnpp_dataset,"anpp_leafnpp_a")
model3[[1]]
model3[[2]]
model3[[3]]
anpp_leafnpp_model <- (lmer(anpp_leafnpp_a~fAPAR_a+vpd_a+PPFD_a+(1|site_a),data=anpp_leafnpp_dataset)) 
r.squaredGLMM(anpp_leafnpp_model)
AIC(anpp_leafnpp_model)
summary(anpp_leafnpp_model)
r.squaredGLMM(anpp_leafnpp_model)
vif_anpp_leafnpp <- vif((lmer(anpp_leafnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset_age)))

#check significance at iteration
summary(lmer(anpp_leafnpp_a~fAPAR_a+(1|site_a),data=anpp_leafnpp_dataset)) 
summary(lmer(anpp_leafnpp_a~fAPAR_a+PPFD_a+(1|site_a),data=anpp_leafnpp_dataset)) 
summary(lmer(anpp_leafnpp_a~fAPAR_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset)) 
summary(lmer(anpp_leafnpp_a~fAPAR_a+PPFD_a+vpd_a+CNrt_a+(1|site_a),data=anpp_leafnpp_dataset)) 
summary(lmer(anpp_leafnpp_a~fAPAR_a+PPFD_a+vpd_a+CNrt_a+Tg_a+(1|site_a),data=anpp_leafnpp_dataset)) 


### BP grassland
#check tnpp grassland
#not filtering any management/non-management! while previous do so
#removing tiandi's grassland
NPP_grassland <- subset(NPP_all,pft=="Grassland" & is.na(Nmin)==TRUE)
grassland_sitemean <- aggregate(NPP_grassland,by=list(NPP_grassland$site), FUN=mean, na.rm=TRUE) 

BP_dataset_grass <- na.omit(grassland_sitemean[,c("tnpp_a","Tg_a","PPFD_a","vpd_a","CNrt_a","fAPAR_a")])
model_g1 <- stepwise_lm(BP_dataset_grass,"tnpp_a")
model_g1[[1]]
model_g1[[2]]
model_g1[[3]]

bp_grass_model <- (lm(tnpp_a~PPFD_a+Tg_a,data=BP_dataset_grass))
summary(bp_grass_model)
r.squaredGLMM(bp_grass_model)

vif_bp_grass <- vif((lm(tnpp_a~Tg_a+PPFD_a+vpd_a+CNrt_a+fAPAR_a,data=BP_dataset_grass)))

#check significance at iteration
summary(lm(tnpp_a~Tg_a,data=BP_dataset_grass))
summary(lm(tnpp_a~PPFD_a+Tg_a,data=BP_dataset_grass))
summary(lm(tnpp_a~PPFD_a+Tg_a+vpd_a,data=BP_dataset_grass))
summary(lm(tnpp_a~PPFD_a+Tg_a+vpd_a+fAPAR_a,data=BP_dataset_grass))
summary(lm(tnpp_a~PPFD_a+Tg_a+vpd_a+fAPAR_a+CNrt_a,data=BP_dataset_grass))

### ANPP grassland
#anpp/tnpp
dim(subset(grassland_sitemean,TNPP_1>0))
dim(subset(grassland_sitemean,TNPP_1>0 & ANPP_2>0))

anpp_tnpp_dataset_grass <- na.omit(grassland_sitemean[,c("anpp_tnpp_a","Tg_a","PPFD_a","vpd_a","CNrt_a","fAPAR_a")])
model_g2 <- stepwise_lm(anpp_tnpp_dataset_grass,"anpp_tnpp_a")
model_g2[[1]]
summary(lm(anpp_tnpp_a~Tg_a,data=anpp_tnpp_dataset_grass))
#non-significant! so alternatively using constant ratio
summary((lm(anpp_a~-1+tnpp_a,data=grassland_sitemean))) # 0.49 for anpp, so 0.51 for bnpp

## leaf Nmass

### Pre-process data
#leaf Nmass
###2. leaf Nmass basing on a site-species model
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
SP_input_sm <- aggregate(SP_input,by=list(SP_input$lon,SP_input$lat), FUN=mean, na.rm=TRUE) #site-mean
summary(SP_input_sm$C_percent)
dim(subset(SP_input_sm,C_percent>0))
#median of value of leaf C = 0.47 or 47%

#remove bahar, as repeated to atkin - and filter a few points with vcmax25<0
SP_input <- subset(SP_input,source!="Bahar et al 2017 New Phytologist")
SP_input <- subset(SP_input,Vcmax25>0)
SP_input$Vcmax.25 <- SP_input$Vcmax25
SP_input$Elevation <- SP_input$z

SP_input2 <- SP_input[,c("lat","lon","Elevation","Vcmax.25","narea","lma")]
sitemean <- aggregate(SP_input2,by=list(SP_input2$lon,SP_input2$lat), FUN=mean, na.rm=TRUE) #site-mean

sitemean$sitename <- paste0("s", 1:nrow(sitemean),sep="") # define sitename (s1,s2..s276)

SP_site1 <- sitemean[,c("lon","lat","sitename")]
SP_final1 <- merge(SP_input,SP_site1,by=c("lat","lon"),all.x=TRUE) #merged sitename to SP data

SP_Vcmax.25 <- aggregate(Vcmax.25~sitename+species,SP_final1,mean) #umol/m2/s
SP_Elevation <- aggregate(Elevation~sitename+species,SP_final1,mean)
SP_narea<- aggregate(narea~sitename+species,SP_final1,mean) # g/m2
SP_lma<- aggregate(lma~sitename+species,SP_final1,mean) # g/m2
SP_lat<- aggregate(lat~sitename+species,SP_final1,mean)
SP_lon<- aggregate(lon~sitename+species,SP_final1,mean)

#merging all observed traits in a site-species dataset.
sitespecies_final <-Reduce(function(x,y) merge(x = x, y = y, by = c("sitename","species"),all.x=TRUE), 
                           list(SP_lon,SP_lat,SP_Elevation,SP_Vcmax.25,
                                SP_narea,SP_lma))

#obtain Nrubisco and Nstructural from this large dataset
#firstly, for site-species data
nmass_a <- sitespecies_final$narea/sitespecies_final$lma
vcmax25_lma_a <- sitespecies_final$Vcmax.25/sitespecies_final$lma
sitename_a <- sitespecies_final$sitename
species_a <- sitespecies_final$species

hist(sitespecies_final$narea) # g/m2
hist(sitespecies_final$lma) # g/m2
hist(sitespecies_final$Vcmax.25) # umol/m2/s

### Fit model
#Fit (Nmass) ~ Ns + Nr * (Vcmax25/LMA) - for site-species data
n1 <- lmer(nmass_a~vcmax25_lma_a + (1|sitename_a)+(1|species_a))
summary(n1)
r.squaredGLMM(n1)

#validation directly

sitemean$pred_nmass <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax.25/sitemean$lma
sitemean$obs_nmass <- sitemean$narea/sitemean$lma

p11 <- analyse_modobs2(sitemean,"pred_nmass","obs_nmass", type = "points",relative=TRUE)$gg + larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest leaf ", N[obs.], " (g g"^-1,")")) +labs(x = ~paste("Forest leaf ", N[pred.], " (g g"^-1,")"))

length(na.omit(sitemean$obs_nmass))

## NRE
### Read and pre-process data
###3. NRE model basing site-mean (lm)
NRE_climate <- read.csv("~/data/NRE_various/NRE_dataset.csv")
NRE_climate$nre_a <- log(NRE_climate$nre/(1-NRE_climate$nre))
NRE_climate$Tg_a <- NRE_climate$Tg
NRE_climate$vpd_a <- log(NRE_climate$vpd)
NRE_climate$PPFD_a <- log(NRE_climate$PPFD)
NRE_climate$ndep_a <- log(NRE_climate$ndep)
NRE_climate2 <- na.omit(NRE_climate[,c("nre_a","vpd_a","Tg_a","PPFD_a","ndep_a")])

### Fit model
stepwise_lm(NRE_climate2,"nre_a")[[1]]
nre_model <- lm(nre_a~Tg_a+vpd_a,data=NRE_climate2)
summary(nre_model)
r.squaredGLMM(nre_model)
length(na.omit(NRE_climate$nre_a))
#validation directly
NRE_climate$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] + summary(nre_model)$coefficients[2,1] *NRE_climate$Tg_a + 
                                      summary(nre_model)$coefficients[3,1] * NRE_climate$vpd_a))))


#visreg

a1 <- ~{
  p1a <- visreg(bp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a2 <- ~{
  p1a <- visreg(bp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln gPPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a3 <- ~{
  p1a <- visreg(bp_model,"CNrt_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln soil C:N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a4 <- ~{
  p1a <- visreg(bp_model,"age_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln age",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a5 <- ~{
  p1a <- visreg(bp_model,"fAPAR_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="fAPAR",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a6 <- ~{
  p1a <- visreg(anpp_tnpp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="logit ABP/BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a7 <- ~{
  p1a <- visreg(anpp_tnpp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="logit ABP/BP",xlab="ln gPPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a8 <- ~{
  p1a <- visreg(anpp_tnpp_model,"CNrt_a",type="contrast")
  plot(p1a,ylab="logit ABP/BP",xlab="ln soil C:N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a9 <- ~{
  p1a <- visreg(anpp_tnpp_model,"age_a",type="contrast")
  plot(p1a,ylab="logit ABP/BP",xlab="ln age",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a10 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"fAPAR_a",type="contrast")
  plot(p1a,ylab="logit leaf-BP/ABP",xlab="fAPAR",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a11 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="logit leaf-BP/ABP",xlab="ln gPPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a12 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"vpd_a",type="contrast")
  plot(p1a,ylab="logit leaf-BP/ABP",xlab="ln D",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a13 <- ~{
  p1a <- visreg(nre_model,"Tg_a",type="contrast")
  plot(p1a,ylab="logit NRE",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a14 <- ~{
  p1a <- visreg(nre_model,"vpd_a",type="contrast")
  plot(p1a,ylab="logit NRE",xlab="ln D",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a15 <- ~{
  p1a <- visreg(bp_grass_model,"Tg_a",type="contrast")
  plot(p1a,ylab="Grassland BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a16 <- ~{
  p1a <- visreg(bp_grass_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="Grassland BP",xlab="ln gPPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

plot_grid(a1,a2,a3,a4,a5,
          a6,a7,a8,a9,white,
          a10,a11,a12,white,white,
          a13,a14,white,white,white,
          a15,a16,white,white,white,
          nrow=5)+white

ggsave(paste("/Users/yunpeng/yunkepeng/CNuptake_MS/output/fig1.jpg",
             sep=""), width = 20, height = 20)

#not run 
#look at validation when using soil C/N
NPP_forest$pred_npp <- summary(bp_model)$coefficients[1,1] +  
  summary(bp_model)$coefficients[2,1] * NPP_forest$Tg_a +
  summary(bp_model)$coefficients[3,1] * NPP_forest$fAPAR_a +
  summary(bp_model)$coefficients[4,1] * NPP_forest$PPFD_a +
  summary(bp_model)$coefficients[5,1] * NPP_forest$CNrt_a+
  summary(bp_model)$coefficients[6,1] * NPP_forest$age_a

NPP_forest$pred_anpp <- NPP_forest$pred_npp * 
  (1/(1 + exp(-(summary(anpp_tnpp_model)$coefficients[1,1]+
                  summary(anpp_tnpp_model)$coefficients[2,1] * NPP_forest$CNrt_a +
                  summary(anpp_tnpp_model)$coefficients[3,1] * NPP_forest$PPFD_a + 
                  summary(anpp_tnpp_model)$coefficients[4,1] * NPP_forest$Tg_a+
                  summary(anpp_tnpp_model)$coefficients[5,1] * NPP_forest$age_a))))

NPP_forest$pred_bnpp <- NPP_forest$pred_npp - NPP_forest$pred_anpp

NPP_forest$pred_lnpp <- NPP_forest$pred_anpp *
  (1/(1 + exp(-(summary(anpp_leafnpp_model)$coefficients[1,1]+
                  summary(anpp_leafnpp_model)$coefficients[2,1]* NPP_forest$fAPAR_a +
                  summary(anpp_leafnpp_model)$coefficients[3,1] * NPP_forest$vpd_a + 
                  summary(anpp_leafnpp_model)$coefficients[4,1] * NPP_forest$PPFD_a))))

NPP_forest$pred_wnpp <- NPP_forest$pred_anpp - NPP_forest$pred_lnpp

NPP_forest$pred_leafnc <- 1/(24.315* NPP_forest$CNrt_a-35.676)

NPP_forest$pred_lnf <- NPP_forest$pred_lnpp*NPP_forest$pred_leafnc
NPP_forest$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                                     summary(nre_model)$coefficients[2,1] *NPP_forest$Tg_a +
                                     summary(nre_model)$coefficients[3,1] * NPP_forest$vpd_a))))

# 100 and 94 are constant, see above for forest wood and root C/N

NPP_forest$pred_wnf <- wood_tissue_percentage*NPP_forest$pred_wnpp/wood_cn_forest
NPP_forest$pred_bnf <- NPP_forest$pred_bnpp/root_cn_forest

NPP_forest$pred_nuptake <- NPP_forest$pred_lnf*(1-NPP_forest$pred_nre)+NPP_forest$pred_wnf+NPP_forest$pred_bnf

NPP_forest$lnf_obs_final <-NPP_forest$NPP.foliage/NPP_forest$CN_leaf_final
NPP_forest$bnf_obs_final  <- NPP_forest$BNPP_1/NPP_forest$CN_root_final
NPP_forest$wnf_obs_final  <- NPP_forest$NPP.wood/NPP_forest$CN_wood_final


p6 <- analyse_modobs2(NPP_forest,"pred_lnf","lnf_obs_final", type = "points",relative=TRUE)$gg+ larger_size+coord_obs_pred()+
  labs(y = ~paste("Forest leaf N ", flux[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest leaf N ", flux[pred.], " (gN m"^-2,"yr"^-1,")"))

p7 <- analyse_modobs2(NPP_forest,"pred_nuptake","Nmin", type = "points",relative=TRUE)$gg+ larger_size+coord_obs_pred()+
  labs(y = ~paste("Net ", minerlization[obs.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")"))

p6
p7
