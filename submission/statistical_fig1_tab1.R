rm(list=ls())
library(lme4)
library(nlme)
library(lmerTest)
library("PerformanceAnalytics")
library(MuMIn)
library(tidyverse)
library(lme4)
library(gridGraphics)
library(visreg)

###1. statistical model of NPP/GPP, ANPP/GPP, leaf NPP/ANPP.
#input data was derived from /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/NPP_statistical_model.Rmd
NPP_statistical <- read.csv("/Users/yunpeng/data/NPP_final/NPP_statistical_forest.csv")

tnpp_gpp_a <- log((NPP_statistical$TNPP_1/NPP_statistical$GPP)/(1-(NPP_statistical$TNPP_1/NPP_statistical$GPP)))
soilCN_a <- log(NPP_statistical$soilCN)
age_a <- log(NPP_statistical$age)
observedfAPAR_a <- NPP_statistical$observedfAPAR
site_a <- NPP_statistical$site

mod_tnpp <- lmer( tnpp_gpp_a ~ soilCN_a + age_a + observedfAPAR_a  + (1|site_a))
summary(mod_tnpp)
r.squaredGLMM(mod_tnpp)
save(mod_tnpp, file = "~/data/NPP_final/statistical_model/mod_tnpp.RData")


anpp_gpp_b <- log((NPP_statistical$ANPP_2/NPP_statistical$GPP)/(1-(NPP_statistical$ANPP_2/NPP_statistical$GPP)))
mod_anpp <- lmer(anpp_gpp_b ~ soilCN_a + age_a + observedfAPAR_a  + (1|site_a))
summary(mod_anpp)
r.squaredGLMM(mod_anpp)
save(mod_anpp, file = "~/data/NPP_final/statistical_model/mod_anpp.RData")


lnpp_data <- subset(NPP_statistical,file=="Sara Vicca"|file=="ForC")
lnpp_anpp_c <- log((lnpp_data$NPP.foliage/lnpp_data$ANPP_2)/(1-(lnpp_data$NPP.foliage/lnpp_data$ANPP_2)))
PPFD_c <- log(lnpp_data$PPFD_gwr)
Tg_c <- lnpp_data$Tg_gwr
vpd_c <- log(lnpp_data$vpd_gwr)
site_c <- lnpp_data$site
mod_lnpp <- lmer(lnpp_anpp_c ~ PPFD_c + Tg_c + vpd_c + (1|site_c))
summary(mod_lnpp)
r.squaredGLMM(mod_lnpp)
save(mod_lnpp, file = "~/data/NPP_final/statistical_model/mod_lnpp.RData")

###2. leaf Nmass basing on a site-species model
SP_input <- read.csv(file="/Users/yunpeng/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
summary(SP_input$C_percent)
#mean of value of leaf C = 0.46 or 46%

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

#Fit (Nmass) ~ Ns + Nr * (Vcmax25/LMA) - for site-species data
n1 <- lmer(nmass_a~vcmax25_lma_a + (1|sitename_a)+(1|species_a))
summary(n1)
r.squaredGLMM(n1)
save(n1, file = "~/data/NPP_final/statistical_model/nmass.RData")

###3. NRE model basing site-mean (lm)
NRE_climate <- read.csv("/Users/yunpeng/data/NPP_final/NRE_statistical_forest.csv")
nre_a <- log(NRE_climate$nre/(1-NRE_climate$nre))
Tg_a <- NRE_climate$Tg
vpd_a <- log(NRE_climate$vpd)
PPFD_a <- log(NRE_climate$PPFD)

nre_model <- lm(nre_a~Tg_a+vpd_a)
summary(nre_model)
save(nre_model, file = "~/data/NPP_final/statistical_model/nre_model_forest.RData")

###4. now forming figure 1 (4 models combing into a panel)
a1 <- ~{
  p1a <- visreg(mod_tnpp,"soilCN_a",type="contrast")
  plot(p1a,ylab="logit (NPP/GPP)",xlab="ln soil C/N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a2 <- ~{
  p1a <- visreg(mod_tnpp,"age_a",type="contrast")
  plot(p1a,ylab=" ",xlab="ln age (years)",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a3 <- ~{
  p1a <- visreg(mod_tnpp,"observedfAPAR_a",type="contrast")
  plot(p1a,ylab=" ",xlab="fAPAR",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a4 <- ~{
  p1a <- visreg(mod_anpp,"soilCN_a",type="contrast")
  plot(p1a,ylab="logit (ANPP/GPP)",xlab="ln soil C/N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a5 <- ~{
  p1a <- visreg(mod_anpp,"age_a",type="contrast")
  plot(p1a,ylab=" ",xlab="ln age (years)",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a6 <- ~{
  p1a <- visreg(mod_anpp,"observedfAPAR_a",type="contrast")
  plot(p1a,ylab=" ",xlab="fAPAR",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a7 <- ~{
  p1a <- visreg(mod_lnpp,"PPFD_c",type="contrast")
  plot(p1a,ylab="logit (leaf-NPP/ANPP)",xlab=(~paste("ln PPFD (", mu, "mol m"^-2,"s"^-1, ")")),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }

a8 <- ~{
  p1a <- visreg(mod_lnpp,"Tg_c",type="contrast")
  plot(p1a,ylab=" ",xlab=(~paste(T[g], " (\u00B0C)")),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  }

a9 <- ~{
  p1a <- visreg(mod_lnpp,"vpd_c",type="contrast")
  plot(p1a,ylab=" ",xlab="D (kPa)",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a10 <- ~{
  p1a <- visreg(nre_model,"Tg_a",type="contrast")
  plot(p1a,ylab="logit (NRE)",xlab=(~paste(T[g], " (\u00B0C)")),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a11 <- ~{
  p1a <- visreg(nre_model,"vpd_a",type="contrast")
  plot(p1a,ylab=" ",xlab="D (kPa)",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}


plot_grid(a1, a2,a3,
          a4, a5,a6,
          a7, a8,a9,
          a10, a11,
          labels = c('(a)','(b)','(c)',
                     '(d)','(e)','(f)',
                     '(g)','(h)', '(i)',
                     '(j)','(k)'),
          nrow=4,label_x = 0.8, label_y = 0.8)

ggsave(paste("~/data/output/fig1.jpg",sep=""), width = 14, height = 16)

#now, prepare table 1
##1. For forst
#NPP/GPP model
summary(mod_tnpp)
#ANPP/GPP model
summary(mod_anpp)
#Leaf-NPP/ANPP model
summary(mod_lnpp)
# Nmass ~ vcmax25/lma model
summary(n1)
# Cmass constant = 46%
SP_input <- read.csv(file="/Users/yunpeng/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
mean(SP_input$C_percent,na.rm=TRUE)
# Root C/N = 94 - as derived from median values of collected samples
NPP_Forest2 <- read.csv("/Users/yunpeng/data/NPP_final/NPP_validation.csv")
summary(NPP_Forest2$CN_root_final) # using median = 94
# Wood C/N = 100 - as derived median values of TRY database
summary(read.csv("~/data/CN_wood/wood_cn.csv")$OrigValueStr)

##2. For Grassland (as orig. data combind from forst_site_org.R and then grassland_simulation.R)
#npp/gpp model = 0.435
NPP_grassland_final4 <- read.csv("/Users/yunpeng/data/NPP_Grassland_final/NPP_grass_validation.csv")
NPP_grassland_final5_gpp_npp_anpp <- aggregate(NPP_grassland_final4,by=list(NPP_grassland_final4$lon,NPP_grassland_final4$lat,NPP_grassland_final4$z), FUN=mean, na.rm=TRUE)
tnpp_grass <- (lm((TNPP_1)~-1+(weightedgpp_all),data=NPP_grassland_final5_gpp_npp_anpp)) #0.435 for using weighted gpp (df = 79)
summary(tnpp_grass)
save(tnpp_grass, file = "~/data/NPP_grassland_final/statistical_model/tnpp_grass.RData")

#anpp/gpp model = 0.228
anpp_grass <- (lm((ANPP_2)~-1+(weightedgpp_all),data=NPP_grassland_final5_gpp_npp_anpp)) #anpp/gpp = 0.228 (df = 289)
summary(anpp_grass)

#leaf c/n model. median = 18
summary(NPP_grassland_final5_gpp_npp_anpp$CN_leaf_final)

#root c/n model median = 41
summary(NPP_grassland_final5_gpp_npp_anpp$CN_root_final)

#NRE = 69%
NRE_Du <- read.csv(file="~/data/NRE_various/NRE_Du/NRE_Du.csv")
NRE_Dong <- read.csv(file="~/data/NRE_various/NRE_Deng/NRE_Deng.csv")

#first - make forest model only
NRE_Du_df <- NRE_Du[,c("lon","lat","NRE","MAT","MAP","VegeType")]
#vegetype =1 is woody ecosystem (all assumed as forest here)
#vegetype = 2 is grassland ecosystem
NRE_Du_df <- subset(NRE_Du_df,VegeType==2)
NRE_Du_df <- aggregate(NRE_Du_df,by=list(NRE_Du_df$lon,NRE_Du_df$lat), FUN=mean, na.rm=TRUE) #site-mean
NRE_Du_df <- NRE_Du_df[,c(3:7)]
head(NRE_Du_df)
dim(NRE_Du_df)

NRE_Dong_df <- NRE_Dong[,c("Longitude","Latitude","NRE.nitrogen.resorption.efficiency.","MAT","MAP","Biome.abbreviation...")]
names(NRE_Dong_df) <- c("lon","lat","NRE","MAT","MAP","biome")
NRE_Dong_df <- subset(NRE_Dong_df,biome=="Grs")
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
summary(NRE_df$NRE)
#median of NRE in grassland is 0.69
