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
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)
library(ggpubr)
#the predicted site values (e.g. predicted npp, anpp..) were already obtained and saved in input csv. 
#But I can show this calculation by statistical model again here
#which is starting from predicted site-level of GPP, Vcmax25 and other predictors (soil C/N, age...) to calculate C-N processes

#Forest NPP validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
NPP_Forest2 <- read.csv("~/data/NPP_final/NPP_validation.csv")
#the predicted site values were already obtained and saved in csv. But I can show again here - starting from predicted site- GPP, Vcmax25 and other predictors
load("~/data/NPP_final/statistical_model/mod_tnpp.RData")
summary(mod_tnpp)
load("~/data/NPP_final/statistical_model/mod_anpp.RData")
summary(mod_anpp)
load("~/data/NPP_final/statistical_model/mod_lnpp.RData")
summary(mod_lnpp)
load("~/data/NPP_final/statistical_model/nmass.RData")
summary(n1)
load("~/data/NPP_final/statistical_model/nre_model_forest.RData")
summary(nre_model)

NPP_Forest2$pred_npp <- NPP_Forest2$pred_gpp_c3 * (1/(1 + exp(-(summary(mod_tnpp)$coefficients[1,1] + 
                                                                summary(mod_tnpp)$coefficients[2,1] * log(NPP_Forest2$CNrt) +
                                                                summary(mod_tnpp)$coefficients[3,1] * log(NPP_Forest2$age) + 
                                                                summary(mod_tnpp)$coefficients[4,1] * NPP_Forest2$fAPAR))))
NPP_Forest2$pred_anpp <- NPP_Forest2$pred_gpp_c3 * (1/(1 + exp(-(summary(mod_anpp)$coefficients[1,1]+
                                                                 summary(mod_anpp)$coefficients[2,1] * log(NPP_Forest2$CNrt) +
                                                                 summary(mod_anpp)$coefficients[3,1] * log(NPP_Forest2$age) + 
                                                                 summary(mod_anpp)$coefficients[4,1] * NPP_Forest2$fAPAR))))
NPP_Forest2$pred_bnpp <- NPP_Forest2$pred_npp - NPP_Forest2$pred_anpp
NPP_Forest2$pred_lnpp <- NPP_Forest2$pred_anpp * (1/(1 + exp(-(summary(mod_lnpp)$coefficients[1,1]+
                                                               summary(mod_lnpp)$coefficients[2,1]* log(NPP_Forest2$PPFD) +
                                                               summary(mod_lnpp)$coefficients[3,1] * (NPP_Forest2$Tg) + 
                                                               summary(mod_lnpp)$coefficients[4,1] * log(NPP_Forest2$vpd)))))
NPP_Forest2$pred_wnpp <- NPP_Forest2$pred_anpp - NPP_Forest2$pred_lnpp
NPP_Forest2$pred_leafnc <- (summary(n1)$coefficients[1,1]/0.46) + (summary(n1)$coefficients[2,1]/0.46) * NPP_Forest2$max_vcmax25/NPP_Forest2$LMA
NPP_Forest2$pred_lnf <- NPP_Forest2$pred_lnpp*NPP_Forest2$pred_leafnc
NPP_Forest2$pred_wnf <- NPP_Forest2$pred_wnpp/100
NPP_Forest2$pred_bnf <- NPP_Forest2$pred_bnpp/94


#Forest leaf Nmass validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
sitemean <- read.csv("~/data/NPP_final/Nmass_validation.csv")
sitemean$pred_leafn <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax25/sitemean$lma
sitemean$obs_leafn <- sitemean$narea/sitemean$lma

#Forest NRE validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
NRE_climate <- read.csv("~/data/NPP_final/NRE_validation.csv")
NRE_climate$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] + summary(nre_model)$coefficients[2,1] *NRE_climate$Tg + summary(nre_model)$coefficients[3,1] * log(NRE_climate$vpd)))))

#Forest Nuptake (Nmin) validation:calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/New_Nuptake_site_simulation.R
Nmin_forest <- read.csv("~/data/NPP_final/Nmin_validation.csv")
Nmin_forest$pred_npp <- Nmin_forest$pred_gpp_c3 * (1/(1 + exp(-(summary(mod_tnpp)$coefficients[1,1] + 
                                                                summary(mod_tnpp)$coefficients[2,1] * log(Nmin_forest$CNrt) +
                                                                summary(mod_tnpp)$coefficients[3,1] * log(Nmin_forest$age) + 
                                                                summary(mod_tnpp)$coefficients[4,1] * Nmin_forest$fAPAR))))
Nmin_forest$pred_anpp <- Nmin_forest$pred_gpp_c3 * (1/(1 + exp(-(summary(mod_anpp)$coefficients[1,1]+
                                                                 summary(mod_anpp)$coefficients[2,1] * log(Nmin_forest$CNrt) +
                                                                 summary(mod_anpp)$coefficients[3,1] * log(Nmin_forest$age) + 
                                                                 summary(mod_anpp)$coefficients[4,1] * Nmin_forest$fAPAR))))
Nmin_forest$pred_bnpp <- Nmin_forest$pred_npp - Nmin_forest$pred_anpp
Nmin_forest$pred_lnpp <- Nmin_forest$pred_anpp * (1/(1 + exp(-(summary(mod_lnpp)$coefficients[1,1]+
                                                               summary(mod_lnpp)$coefficients[2,1]* log(Nmin_forest$PPFD) +
                                                               summary(mod_lnpp)$coefficients[3,1] * (Nmin_forest$Tg) + 
                                                               summary(mod_lnpp)$coefficients[4,1] * log(Nmin_forest$vpd)))))
Nmin_forest$pred_wnpp <- Nmin_forest$pred_anpp - Nmin_forest$pred_lnpp
Nmin_forest$pred_leafnc <- (summary(n1)$coefficients[1,1]/0.46) + (summary(n1)$coefficients[2,1]/0.46) * Nmin_forest$max_vcmax25/Nmin_forest$LMA
Nmin_forest$pred_lnf <- Nmin_forest$pred_lnpp*Nmin_forest$pred_leafnc
Nmin_forest$pred_wnf <- Nmin_forest$pred_wnpp/100
Nmin_forest$pred_bnf <- Nmin_forest$pred_bnpp/94
Nmin_forest$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] + summary(nre_model)$coefficients[2,1] *Nmin_forest$Tg + summary(nre_model)$coefficients[3,1] * log(Nmin_forest$vpd)))))
Nmin_forest$pred_nuptake <- Nmin_forest$pred_lnf*(1-Nmin_forest$pred_nre) + Nmin_forest$pred_wnf +Nmin_forest$pred_bnf  

#conduct additional analysis
#the relationship between c/n of leaf, root and wood was available in (because we newly need Schulze's tissues data for wood C/N!)
#load(file = "/Users/yunpeng/data/NPP_final/Forest_site_simulation.Rdata")
#summary(lm(NPP_Forest$CN_root_final~NPP_Forest$CN_leaf_org))
#summary(lm(NPP_Forest$CN_wood_final~NPP_Forest$CN_leaf_org))

#conduct bias analysis between BP[obs] - BP[pred] vs. soil C/N (expect a negative relationship - see Sara's email!)
My_Theme = theme(
  axis.title.x = element_text(size = 15),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 15),
  axis.text.y = element_text(size = 20))+theme_classic()

NPP_Forest2$bias_bp <- (NPP_Forest2$TNPP_1 - NPP_Forest2$pred_npp)
summary(lm(NPP_Forest2$bias_bp~NPP_Forest2$soilCN))
summary(lmer(bias_bp~soilCN+(1|site),data=NPP_Forest2))
NPP_Forest2$bias_bpe <- NPP_Forest2$TNPP_1/NPP_Forest2$GPP - NPP_Forest2$pred_npp/NPP_Forest2$pred_gpp_c3
NPP_Forest2$bias_gpp <- NPP_Forest2$GPP - NPP_Forest2$pred_gpp_c3

q1 <- ggplot(data=NPP_Forest2, aes(x=soilCN, y=bias_bpe)) +
  geom_point()+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Bias of BPE")) +labs(x = ~paste("soil C/N "))+stat_cor(method = "pearson")

q2 <- ggplot(data=NPP_Forest2, aes(x=soilCN, y=bias_bp)) +
  geom_point()+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Bias of BP")) +labs(x = ~paste("soil C/N "))+stat_cor(method = "pearson")

q3 <- ggplot(data=NPP_Forest2, aes(x=soilCN, y=bias_gpp)) +
  geom_point()+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Bias of GPP")) +labs(x = ~paste("soil C/N "))+stat_cor(method = "pearson")

plot_grid(q1,q2,q3)

ggsave(paste("~/data/output/si.jpg",sep=""),width = 15, height = 16)

#Grassland NPP and leaf N validation validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Grassland_simulation.R
load(file = "~/data/NPP_Grassland_final/statistical_model/tnpp_grass.RData")
mod_tnpp_grass<- tnpp_grass
summary(mod_tnpp_grass)
load(file = "~/data/NPP_Grassland_final/statistical_model/anpp_grass.RData")
mod_anpp_grass <- anpp_grass
summary(mod_anpp_grass)

NPP_grassland_final4 <- read.csv("~/data/NPP_Grassland_final/NPP_grass_validation.csv")
NPP_grassland_final4$pred_npp <- NPP_grassland_final4$weightedgpp_all * summary(mod_tnpp_grass)$coef[1,1]
NPP_grassland_final4$pred_anpp <- NPP_grassland_final4$weightedgpp_all * summary(mod_anpp_grass)$coef[1,1]
NPP_grassland_final4$pred_lnf <- NPP_grassland_final4$pred_anpp/18
NPP_grassland_final4$pred_bnpp <- NPP_grassland_final4$pred_npp - NPP_grassland_final4$pred_anpp
NPP_grassland_final4$pred_bnf <- NPP_grassland_final4$pred_bnpp/41

#validation results for figure 2 and table s2 (plots, R2, RMSE in percentage)
p1 <- ggplot(data=NPP_Forest2, aes(x=pred_gpp_c3, y=GPP)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest ", GPP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", GPP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,5000)+ylim(0,5000)
a1 <- (lm(GPP~pred_gpp_c3,NPP_Forest2))
summary(a1)
mean_value <- mean(subset(NPP_Forest2,pred_gpp_c3>0 & GPP>0)$GPP)
sqrt(mean(a1$residuals^2))/mean_value
sqrt(mean(a1$residuals^2))

p2 <- ggplot(data=NPP_Forest2, aes(x=pred_npp, y=TNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest ", BP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,2000)+ylim(0,2000)

a2 <- (lm(TNPP_1~pred_npp,NPP_Forest2))
summary(a2)
mean_value <- mean(subset(NPP_Forest2,pred_npp>0 & TNPP_1>0)$TNPP_1)
sqrt(mean(a2$residuals^2))/mean_value
sqrt(mean(a2$residuals^2))

p3 <-ggplot(data=NPP_Forest2, aes(x=pred_anpp, y=ANPP_2)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest ", ANPP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", ANPP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,2000)+ylim(0,2000)
a3 <- (lm(ANPP_2~pred_anpp,NPP_Forest2))
summary(a3)
mean_value <- mean(subset(NPP_Forest2,pred_anpp>0 & ANPP_2>0)$ANPP_2)
sqrt(mean(a3$residuals^2))/mean_value
sqrt(mean(a3$residuals^2))

p4 <-ggplot(data=NPP_Forest2, aes(x=pred_lnpp, y=NPP.foliage)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest leaf ", NPP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest leaf ", NPP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,1000)+ylim(0,1000)

a4 <- (lm(NPP.foliage~pred_lnpp,NPP_Forest2))
summary(a4)
mean_value <- mean(subset(NPP_Forest2,pred_lnpp>0 & NPP.foliage>0)$NPP.foliage)
sqrt(mean(a4$residuals^2))/mean_value
sqrt(mean(a4$residuals^2))

p5 <-ggplot(data=NPP_Forest2, aes(x=pred_wnpp, y=NPP.wood)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest wood ", NPP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest wood ", NPP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,1500)+ylim(0,1500)

a5 <- (lm(NPP.wood~pred_wnpp,NPP_Forest2))
summary(a5)
mean_value <- mean(subset(NPP_Forest2,pred_wnpp>0 & NPP.wood>0)$NPP.wood)
sqrt(mean(a5$residuals^2))/mean_value
sqrt(mean(a5$residuals^2))

p6 <-ggplot(data=NPP_Forest2, aes(x=pred_bnpp, y=BNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest ", BNPP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest ", BNPP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,1000)+ylim(0,1000)
a6 <- (lm(BNPP_1~pred_bnpp,NPP_Forest2))
summary(a6)
mean_value <- mean(subset(NPP_Forest2,pred_bnpp>0 & BNPP_1>0)$BNPP_1)
sqrt(mean(a6$residuals^2))/mean_value
sqrt(mean(a6$residuals^2))

p7 <-ggplot(data=NPP_Forest2, aes(x=pred_lnf, y=lnf_obs_org)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest leaf N ", flux[pred.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest leaf N ", flux[obs.], " (gN m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,35)+ylim(0,35)
a7 <- (lm(lnf_obs_org~pred_lnf,NPP_Forest2))
summary(a7)
mean_value <- mean(subset(NPP_Forest2,pred_lnf>0 & lnf_obs_org>0)$lnf_obs_org)
sqrt(mean(a7$residuals^2))/mean_value
sqrt(mean(a7$residuals^2))

p8 <- ggplot(data=sitemean, aes(x=pred_leafn, y=obs_leafn)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest leaf ", N[pred.], " (g g"^-1,")")) +labs(x = ~paste("Forest leaf ", N[obs.], " (g g"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+xlim(0,0.06)+ylim(0,0.06)
a8 <- (lm(obs_leafn~pred_leafn,sitemean))
summary(a8)
mean_value <- mean(subset(sitemean,pred_leafn>0 & obs_leafn>0)$obs_leafn)
sqrt(mean(a8$residuals^2))/mean_value
sqrt(mean(a8$residuals^2))

p9 <- ggplot(data=NRE_climate, aes(x=pred_nre, y=nre)) + #xlim(c(0.25,1))+ylim(c(0.25,1))+
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest ", NRE[pred.])) +labs(x = ~paste("Forest ", NRE[obs.]))+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+xlim(0,1)+ylim(0,1)
a9 <- (lm(nre~pred_nre,NRE_climate))
summary(a9)
mean_value <- mean(subset(NRE_climate,pred_nre>0 & nre>0)$nre)
sqrt(mean(a9$residuals^2))/mean_value
sqrt(mean(a9$residuals^2))
#take care about Nmin - might be cancelling range of Nmin?
#there are 9 samples (within 3 plots) having Nmin >25. But these points are not be able to get to pred Nmin - perhaps because coord are wrong!
p10 <- ggplot(data=Nmin_forest, aes(x=pred_nuptake, y=Nmin)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Forest N ", uptake[pred.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Forest N ", uptake[obs.], " (gN m"^-2,"yr"^-1,")"))+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+xlim(0,25)+ylim(0,25)+theme(  panel.border = element_rect(colour = "black", fill=NA, size=5))
a10 <- (lm(Nmin~pred_nuptake,Nmin_forest))
summary(a10)
mean_value <- mean(subset(Nmin_forest,pred_nuptake>0 & Nmin>0)$Nmin)
sqrt(mean(a10$residuals^2))/mean_value
sqrt(mean(a10$residuals^2))

#plot points of N minerlization rate as raised by Keith.
My_Theme2 = theme(
  axis.title.x = element_text(size = 15),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 15),
  axis.text.y = element_text(size = 20))+ 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

Nmin_forest
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(Nmin_forest$lon,Nmin_forest$lat, col="red", pch=16,cex=1)

p11 <- ggplot(NPP_grassland_final4, aes(x=weightedgpp_all, y=GPP)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Grassland ", GPP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland ", GPP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,3000)+ylim(0,3000)
a11 <- (lm(GPP~weightedgpp_all,NPP_grassland_final4))
summary(a11)
mean_value <- mean(subset(NPP_grassland_final4,weightedgpp_all>0 & GPP>0)$GPP)
sqrt(mean(a11$residuals^2))/mean_value
sqrt(mean(a11$residuals^2))

p12 <-ggplot(NPP_grassland_final4, aes(x=pred_npp, y=TNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Grassland ", BP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland ", BP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,2500)+ylim(0,2500)
a12 <- (lm(TNPP_1~pred_npp,NPP_grassland_final4))
summary(a12)
mean_value <- mean(subset(NPP_grassland_final4,pred_npp>0 & TNPP_1>0)$TNPP_1)
sqrt(mean(a12$residuals^2))/mean_value
sqrt(mean(a12$residuals^2))

p13 <-ggplot(NPP_grassland_final4, aes(x=pred_anpp, y=ANPP_2)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  My_Theme+labs(y = ~paste("Grassland ", ANPP[pred.], " (gC m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland ", ANPP[obs.], " (gC m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,2000)+ylim(0,2000)
a13 <- (lm(ANPP_2~pred_anpp,NPP_grassland_final4))
summary(a13)
mean_value <- mean(subset(NPP_grassland_final4,pred_anpp>0 & ANPP_2>0)$ANPP_2)
sqrt(mean(a13$residuals^2))/mean_value
sqrt(mean(a13$residuals^2))

p14 <-ggplot(NPP_grassland_final4, aes(x=pred_lnf, y=lnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+ #xlim(c(0,10))+
  My_Theme+labs(y = ~paste("Grassland leaf N ", flux[pred.], " (gN m"^-2,"yr"^-1,")")) +labs(x = ~paste("Grassland leaf N ", flux[obs.], " (gN m"^-2,"yr"^-1,")")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,40)+ylim(0,40)
a14 <- (lm(lnf_obs_final~pred_lnf,NPP_grassland_final4))
summary(a14)
mean_value <- mean(subset(NPP_grassland_final4,pred_lnf>0 & lnf_obs_final>0)$lnf_obs_final)
sqrt(mean(a14$residuals^2))/mean_value
sqrt(mean(a14$residuals^2))

plot_grid(p1,p2,p3,
          p11,p12,p13,
          p4,p5,p6,
          p7,p8,p9,
          p14,p10, labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)'),
          nrow=5,ncol=3,
          label_x = 0.9,label_y=0.92)

ggsave(paste("~/data/output/fig2.jpg",sep=""),width = 15, height = 20)

NPP_Forest2$BPE_obs <- NPP_Forest2$TNPP_1/NPP_Forest2$GPP
NPP_Forest2$BPE_pred <- NPP_Forest2$pred_npp/NPP_Forest2$pred_gpp_c3

ggplot(data=NPP_Forest2, aes(x=BPE_pred, y=BPE_obs,color=lat)) +geom_point()+
  scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", space ="Lab" )
ggplot(data=NPP_Forest2, aes(x=lat, y=BPE_obs)) +geom_point()+geom_smooth()
ggplot(data=NPP_Forest2, aes(x=lat, y=BPE_pred)) +geom_point()+geom_smooth()

  
ggplot(data=NPP_Forest2, aes(x=pred_npp, y=TNPP_1)) +geom_point()+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)
ggplot(data=NPP_Forest2, aes(x=pred_gpp_c3, y=GPP)) +geom_point()+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)

#additional soil P collection analysis suggested by sara vicca, metteo and beni
#dim(NPP_Forest2)
#NPP_Forest2_sm <- aggregate(NPP_Forest2,by=list(NPP_Forest2$lon,NPP_Forest2$lat), mean,na.rm=TRUE)
#NPP_input <- NPP_Forest2_sm[,c("lon","lat")]
#NPP_input$sitename <- paste("a",1:nrow(NPP_input),sep="")
#now, doing soil P collection in desktop (takes time commented out here)
#devtools::load_all("/Users/yunpeng/yunkepeng/gcme/pmodel/ingestr/")
#df_layers <- tibble(layer = 1:8, bottom = c(4.5, 9.1, 16.6, 28.9, 49.3, 82.9, 138.3, 229.6)) %>% mutate(top = lag(bottom)) %>% mutate(top = ifelse(is.na(top), 0, top))
#settings_gsde <- list(varnam = c("PBR"), layer = 1:3)
#NPP_input$soil_P <- NA
#for (i in 1:nrow(NPP_input)){
#  df_gsde <- ingest_bysite(sitename  = NPP_input$sitename[i],source = "gsde",lon = NPP_input$lon[i],lat= NPP_input$lat[i],settings  = settings_gsde,dir = "/Volumes/My Passport/data/soil/shangguan")
#  value <- as.numeric(as.data.frame(df_gsde$data[[1]]))
#  NPP_input[i,"soil_P"] <- value
#  print(i)}
#csvfile <- paste("~/data/NPP_final/soilP.csv")
#write.csv(NPP_input, csvfile, row.names = TRUE)
soilP <- read.csv("~/data/NPP_final/soilP.csv")

NPP_Forest2_soilP <- merge(NPP_Forest2,soilP,by=c("lon","lat"),all.x=TRUE)
NPP_Forest2_soilP$bias_bpe <- NPP_Forest2_soilP$TNPP_1/NPP_Forest2_soilP$GPP - NPP_Forest2_soilP$pred_npp/NPP_Forest2_soilP$pred_gpp_c3

ggplot(data=NPP_Forest2_soilP, aes(x=soil_P, y=bias_bpe)) +geom_point()+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)

#add a few measured soil P data?
#add site-level soil C/N
Sara_CN <- read.csv(file="/Users/yunpeng/data/NPP_Yunke/NPP_Vicca/orig/References_Yunke_soilCN.csv")
Sara_CN <- Sara_CN[,c("Plot.name","Total.P..mg.P.kg.1.")]
Sara_CN$soilP <- (as.numeric(gsub(",",".",Sara_CN[,2])))
hist(Sara_CN$soilP)
Sara_CN_site <- aggregate(Sara_CN,by=list(Sara_CN$Plot.name), FUN=mean, na.rm=TRUE)
Sara_CN_site <- subset(Sara_CN_site,soilP>0)
Sara_CN_site$Plot <- Sara_CN_site$Group.1
Sara_CN_site <- Sara_CN_site[,c("Plot","soilP")]
names(Sara_CN_site) <- c("site","soilP_measured")
NPP_Forest2_soilP_measured <- merge(NPP_Forest2_soilP,Sara_CN_site,by=c("site"),all.x=TRUE)
summary(NPP_Forest2_soilP_measured)
dim(NPP_Forest2_soilP_measured) #21 sites available
plot(NPP_Forest2_soilP_measured$soilP_measured~NPP_Forest2_soilP_measured$soil_P)

ggplot(data=NPP_Forest2_soilP_measured, aes(x=soilP_measured, y=bias_bpe)) +geom_point()+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)
summary(lm(NPP_Forest2_soilP_measured$bias_bpe~NPP_Forest2_soilP_measured$soilP_measured))

NPP_statistical <- read.csv("~/data/NPP_final/NPP_statistical_forest.csv")
NPP_statistical <- merge(NPP_statistical,Sara_CN_site,by=c("site"),all.x=TRUE)
tnpp_gpp_a <- log((NPP_statistical$TNPP_1/NPP_statistical$GPP)/(1-(NPP_statistical$TNPP_1/NPP_statistical$GPP)))
soilCN_a <- log(NPP_statistical$soilCN)
age_a <- log(NPP_statistical$age)
observedfAPAR_a <- NPP_statistical$observedfAPAR
site_a <- NPP_statistical$site
P_a <- log(NPP_statistical$soilP_measured)

geom_hex()

mod_tnpp <- lmer( tnpp_gpp_a ~ soilCN_a + age_a + observedfAPAR_a + (1|site_a))
summary(mod_tnpp)
mod_tnpp <- lmer( tnpp_gpp_a ~ P_a  + (1|site_a))
summary(mod_tnpp)
