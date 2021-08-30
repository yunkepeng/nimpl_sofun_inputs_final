rm(list=ls())
#Forest NPP validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
NPP_Forest2 <- read.csv("/Users/yunpeng/data/NPP_final/NPP_validation.csv")

#Forest leaf Nmass validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
sitemean <- read.csv("/Users/yunpeng/data/NPP_final/Nmass_validation.csv")

#Forest NRE validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
NRE_climate <- read.csv("/Users/yunpeng/data/NPP_final/NRE_validation.csv")

#Forest Nuptake (Nmin) validation:calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/New_Nuptake_site_simulation.R
Nmin_forest <- read.csv("/Users/yunpeng/data/NPP_final/Nmin_validation.csv")

#Grassland NPP and leaf N validation validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Grassland_simulation.R
NPP_grassland_final4 <- read.csv("/Users/yunpeng/data/NPP_Grassland_final/NPP_grass_validation.csv")

My_Theme = theme(
  axis.title.x = element_text(size = 15),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 15),
  axis.text.y = element_text(size = 20))+theme_classic()

#validation results for figure 2 and table s2 (plots, R2, RMSE in percentage)
p1 <- ggplot(data=NPP_Forest2, aes(x=pred_gpp_c3, y=GPP)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest ", GPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest ", GPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a1 <- (lm(GPP~pred_gpp_c3,NPP_Forest2))
summary(a1)
mean_value <- mean(subset(NPP_Forest2,pred_gpp_c3>0 & GPP>0)$GPP)
sqrt(mean(a1$residuals^2))/mean_value

p2 <- ggplot(data=NPP_Forest2, aes(x=pred_npp, y=TNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest ", NPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest ", NPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a2 <- (lm(TNPP_1~pred_npp,NPP_Forest2))
summary(a2)
mean_value <- mean(subset(NPP_Forest2,pred_npp>0 & TNPP_1>0)$TNPP_1)
sqrt(mean(a2$residuals^2))/mean_value

p3 <-ggplot(data=NPP_Forest2, aes(x=pred_anpp, y=ANPP_2)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest ", ANPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest ", ANPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a3 <- (lm(ANPP_2~pred_anpp,NPP_Forest2))
summary(a3)
mean_value <- mean(subset(NPP_Forest2,pred_anpp>0 & ANPP_2>0)$ANPP_2)
sqrt(mean(a3$residuals^2))/mean_value

p4 <-ggplot(data=NPP_Forest2, aes(x=pred_lnpp, y=NPP.foliage)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest leaf ", NPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest leaf ", NPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a4 <- (lm(NPP.foliage~pred_lnpp,NPP_Forest2))
summary(a4)
mean_value <- mean(subset(NPP_Forest2,pred_lnpp>0 & NPP.foliage>0)$NPP.foliage)
sqrt(mean(a4$residuals^2))/mean_value

p5 <-ggplot(data=NPP_Forest2, aes(x=pred_wnpp, y=NPP.wood)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest wood ", NPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest wood ", NPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a5 <- (lm(NPP.wood~pred_wnpp,NPP_Forest2))
summary(a5)
mean_value <- mean(subset(NPP_Forest2,pred_wnpp>0 & NPP.wood>0)$NPP.wood)
sqrt(mean(a5$residuals^2))/mean_value

p6 <-ggplot(data=NPP_Forest2, aes(x=pred_bnpp, y=BNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest ", BNPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest ", BNPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a6 <- (lm(BNPP_1~pred_bnpp,NPP_Forest2))
summary(a6)
mean_value <- mean(subset(NPP_Forest2,pred_bnpp>0 & BNPP_1>0)$BNPP_1)
sqrt(mean(a6$residuals^2))/mean_value

p7 <-ggplot(data=NPP_Forest2, aes(x=pred_lnf, y=lnf_obs_org)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest leaf N ", flux[pred.], " (gN m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest leaf N ", flux[obs.], " (gN m"^-2,"s"^-1,")")) 
a7 <- (lm(lnf_obs_org~pred_lnf,NPP_Forest2))
summary(a7)
mean_value <- mean(subset(NPP_Forest2,pred_lnf>0 & lnf_obs_org>0)$lnf_obs_org)
sqrt(mean(a7$residuals^2))/mean_value

p8 <- ggplot(data=sitemean, aes(x=pred_leafn, y=obs_leafn)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest leaf ", N[pred.], " (g m"^-2,")")) +labs(x = ~paste("Forest leaf ", N[obs.], " (g m"^-2,")")) 
a8 <- (lm(obs_leafn~pred_leafn,sitemean))
summary(a8)
mean_value <- mean(subset(sitemean,pred_leafn>0 & obs_leafn>0)$obs_leafn)
sqrt(mean(a8$residuals^2))/mean_value

p9 <- ggplot(data=NRE_climate, aes(x=pred_nre, y=nre)) + xlim(c(0.25,1))+ylim(c(0.25,1))+
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Forest ", NRE[pred.])) +labs(x = ~paste("Forest ", NRE[obs.]))
a9 <- (lm(nre~pred_nre,NRE_climate))
summary(a9)
mean_value <- mean(subset(NRE_climate,pred_nre>0 & nre>0)$nre)
sqrt(mean(a9$residuals^2))/mean_value

p10 <- ggplot(data=Nmin_forest, aes(x=pred_nuptake, y=Nmin)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+xlim(0,25)+ylim(0,25)+
  My_Theme+labs(y = ~paste("Forest N ", uptake[pred.], " (gN m"^-2,"s"^-1,")")) +labs(x = ~paste("Forest N ", uptake[obs.], " (gN m"^-2,"s"^-1,")"))
a10 <- (lm(Nmin~pred_nuptake,Nmin_forest))
summary(a10)
mean_value <- mean(subset(Nmin_forest,pred_nuptake>0 & Nmin>0)$Nmin)
sqrt(mean(a10$residuals^2))/mean_value

p11 <- ggplot(NPP_grassland_final4, aes(x=weightedgpp_all, y=GPP)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Grassland ", GPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Grassland ", GPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a11 <- (lm(GPP~weightedgpp_all,NPP_grassland_final4))
summary(a11)
mean_value <- mean(subset(NPP_grassland_final4,weightedgpp_all>0 & GPP>0)$GPP)
sqrt(mean(a11$residuals^2))/mean_value

p12 <-ggplot(NPP_grassland_final4, aes(x=pred_npp, y=TNPP_1)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Grassland ", NPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Grassland ", NPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a12 <- (lm(TNPP_1~pred_npp,NPP_grassland_final4))
summary(a12)
mean_value <- mean(subset(NPP_grassland_final4,pred_npp>0 & TNPP_1>0)$TNPP_1)
sqrt(mean(a12$residuals^2))/mean_value


p13 <-ggplot(NPP_grassland_final4, aes(x=pred_anpp, y=ANPP_2)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+
  My_Theme+labs(y = ~paste("Grassland ", ANPP[pred.], " (gC m"^-2,"s"^-1,")")) +labs(x = ~paste("Grassland ", ANPP[obs.], " (gC m"^-2,"s"^-1,")")) 
a13 <- (lm(ANPP_2~pred_anpp,NPP_grassland_final4))
summary(a13)
mean_value <- mean(subset(NPP_grassland_final4,pred_anpp>0 & ANPP_2>0)$ANPP_2)
sqrt(mean(a13$residuals^2))/mean_value

p14 <-ggplot(NPP_grassland_final4, aes(x=pred_lnf, y=lnf_obs_final)) +
  geom_point()+geom_abline(intercept=0,slope=1)+geom_smooth(method = "lm", se = TRUE)+ xlim(c(0,10))+
  My_Theme+labs(y = ~paste("Grassland leaf N ", flux[pred.], " (gN m"^-2,"s"^-1,")")) +labs(x = ~paste("Grassland leaf N ", flux[obs.], " (gN m"^-2,"s"^-1,")")) 
a14 <- (lm(lnf_obs_final~pred_lnf,NPP_grassland_final4))
summary(a14)
mean_value <- mean(subset(NPP_grassland_final4,pred_lnf>0 & lnf_obs_final>0)$lnf_obs_final)
sqrt(mean(a14$residuals^2))/mean_value

plot_grid(p1, p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)',
                                                                      '(i)','(j)','(k)','(l)','(m)','(n)'),label_x = 0.9)

ggsave(paste("~/data/output/fig2.jpg",sep=""),width = 15, height = 16)
