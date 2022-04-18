#Forest NPP validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
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
library(caret)
library(recipes)
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")
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

rm(list=ls())

#1. input complete dataset and do some re-processing
stepwise <- function(df_input,target_var){
  #-----------------------------------------------------------------------
  # Input:  whole dataframe and target variable
  #assume that site_a is the only random factor
  #-----------------------------------------------------------------------
  target <- target_var
  df <- df_input
  
  preds <- df %>% dplyr::select(-c(target,site_a)) %>% 
    names()
  
  r_list <- c()
  #For loop functions, include all predictor's r2 at the end
  for (var in preds){
    forml <- paste( 'lmer(', target, '~', var, '+(1|site_a), data = df)')
    fit_lin <- eval(parse(text = forml)) 
    rsq <- r.squaredGLMM(fit_lin)[1]
    r_list <- c(r_list,rsq)
  }
  
  #convert to a dataframe, including all r2
  All_rsquare <- data.frame (
    preds = factor(preds,levels=preds), 
    rsq = r_list)
  
  #select max r2 in all predictors
  max(r_list)
  
  new_All_rsquare <- All_rsquare %>% 
    # desc orders from largest to smallest
    arrange(desc(rsq))
  
  #2. stepwise regression selection
  
  ## list
  list_aic <- list()
  list_bic <- list()
  list_R <- list()
  list_variable <- list()
  
  # predictors retained in the model firstly
  preds_retained <- as.character(new_All_rsquare[1,1])
  preds_candidate <- preds[-which(preds == preds_retained)] 
  
  
  for (a in 1:(length(preds)-1)){
    rsq_candidates <- c()
    linmod_candidates <- list()
    for (i in 1:length(preds_candidate)){
      pred_add <- c(preds_retained, preds_candidate[i])
      forml  <- paste( 'lmer(', target, '~', paste(pred_add, collapse = '+'), '+(1|site_a), data = df)')
      # create a function and make its format available to output in for loop
      fit_lin <- eval(parse(text = forml))
      linmod_candidates[[ i ]] <- fit_lin
      # obtain multiple r2 at each selection, and find the best one at the end
      rsq <- r.squaredGLMM(fit_lin)[1]
      rsq_candidates[i] <- rsq
    }
    pred_max <- preds_candidate[ which.max(rsq_candidates) ]
    # include best factors in retained factor
    preds_retained <- c(preds_retained, pred_max)
    list_variable[[a]] <- pred_max 
    # include AIC, BIC, adjusted R2, R2, cross-validated R2 and RMSE at each k 
    list_aic[[  a ]] <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = df)'))))
    
    list_bic[[ a ]] <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = df)'))))
    
    list_R[[ a ]] <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = df)'))))[1]
    preds_candidate <- preds_candidate[-which(preds_candidate == pred_max)]
  }
  
  
  R_null <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = df)'))))[1]
  AIC_null <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = df)'))))
  BIC_null <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = df)'))))
  variable_null <- preds_retained[1]
  
  R_all <- round(as.numeric(c(R_null,list_R)),2)
  AIC_all <- round(as.numeric(c(AIC_null,list_aic)),2)
  BIC_all <- round(as.numeric(c(BIC_null,list_bic)),2)
  variable_all <- (as.character(c(variable_null,list_variable)))
  
  df1 <- as.data.frame(cbind(variable_all,R_all,AIC_all,BIC_all))
  
  #Adjusted-R
  p1 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = R_all)) 
  #AIC
  p2 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = AIC_all)) 
  #BIC
  p3 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = BIC_all))
  
  output_list <- list(p1,p2,p3)
  
  return(output_list)
  #-----------------------------------------------------------------------
  # Output: four figures 
  #-----------------------------------------------------------------------
}

#stepwise lm
stepwise_lm <- function(df_input,target_var){
  #-----------------------------------------------------------------------
  # Input:  whole dataframe and target variable
  #-----------------------------------------------------------------------
  target <- target_var
  df <- df_input
  
  preds <- df %>% dplyr::select(-c(target)) %>% 
    names()
  
  r_list <- c()
  #For loop functions, include all predictor's r2 at the end
  for (var in preds){
    forml <- paste( 'lm(', target, '~', var, ', data = df)')
    fit_lin <- eval(parse(text = forml)) 
    rsq <- r.squaredGLMM(fit_lin)[1]
    r_list <- c(r_list,rsq)
  }
  
  #convert to a dataframe, including all r2
  All_rsquare <- data.frame (
    preds = factor(preds,levels=preds), 
    rsq = r_list)
  
  #select max r2 in all predictors
  max(r_list)
  
  new_All_rsquare <- All_rsquare %>% 
    # desc orders from largest to smallest
    arrange(desc(rsq))
  
  #2. stepwise regression selection
  
  ## list
  list_aic <- list()
  list_bic <- list()
  list_R <- list()
  list_variable <- list()
  
  # predictors retained in the model firstly
  preds_retained <- as.character(new_All_rsquare[1,1])
  preds_candidate <- preds[-which(preds == preds_retained)] 
  
  
  for (a in 1:(length(preds)-1)){
    rsq_candidates <- c()
    linmod_candidates <- list()
    for (i in 1:length(preds_candidate)){
      pred_add <- c(preds_retained, preds_candidate[i])
      forml  <- paste( 'lm(', target, '~', paste(pred_add, collapse = '+'), ', data = df)')
      # create a function and make its format available to output in for loop
      fit_lin <- eval(parse(text = forml))
      linmod_candidates[[ i ]] <- fit_lin
      # obtain multiple r2 at each selection, and find the best one at the end
      rsq <- r.squaredGLMM(fit_lin)[1]
      rsq_candidates[i] <- rsq
    }
    pred_max <- preds_candidate[ which.max(rsq_candidates) ]
    # include best factors in retained factor
    preds_retained <- c(preds_retained, pred_max)
    list_variable[[a]] <- pred_max 
    # include AIC, BIC, adjusted R2, R2, cross-validated R2 and RMSE at each k 
    list_aic[[  a ]] <- AIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained, collapse = '+'),  ', data = df)'))))
    
    list_bic[[ a ]] <- BIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained, collapse = '+'),  ', data = df)'))))
    
    list_R[[ a ]] <- r.squaredGLMM(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained, collapse = '+'),  ', data = df)'))))[1]
    preds_candidate <- preds_candidate[-which(preds_candidate == pred_max)]
  }
  
  
  R_null <- r.squaredGLMM(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained[1], collapse = '+'),  ', data = df)'))))[1]
  AIC_null <- AIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained[1], collapse = '+'),  ', data = df)'))))
  BIC_null <- BIC(eval(parse(text = paste( 'lm(', target, '~', paste(preds_retained[1], collapse = '+'),  ', data = df)'))))
  variable_null <- preds_retained[1]
  
  R_all <- round(as.numeric(c(R_null,list_R)),2)
  AIC_all <- round(as.numeric(c(AIC_null,list_aic)),2)
  BIC_all <- round(as.numeric(c(BIC_null,list_bic)),2)
  variable_all <- (as.character(c(variable_null,list_variable)))
  
  df1 <- as.data.frame(cbind(variable_all,R_all,AIC_all,BIC_all))
  
  #Adjusted-R
  p1 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = R_all)) 
  #AIC
  p2 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = AIC_all)) 
  #BIC
  p3 <- ggplot() + 
    geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = BIC_all))
  
  output_list <- list(p1,p2,p3)
  
  return(output_list)
  #-----------------------------------------------------------------------
  # Output: four figures 
  #-----------------------------------------------------------------------
}

#remove alpha throughout all analyses because 
#(1) D and alpha shows contridictory results for NUE 
#(2) D and alpha repeated in model selection

NPP_all <- read.csv("~/data/NPP_Yunke/NPP_Nmin_dataset_with_predictors.csv")

#remove tiandi's grassland since their bnpp and bp were considered to be not used
#tiandi's bnpp and bp were alredy removed in Forest_site_orig.R
NPP_all$TNPP_1[NPP_all$file=="Tiandi Grassland"] <- NA
NPP_all$ANPP_2[NPP_all$file=="Tiandi Grassland"] <- NA
NPP_all$BNPP_1[NPP_all$file=="Tiandi Grassland"] <- NA
NPP_all$NPP.foliage[NPP_all$file=="Tiandi Grassland"] <- NA
NPP_all$lnf_obs_final[NPP_all$file=="Tiandi Grassland"] <- NA
NPP_all$bnf_obs_final[NPP_all$file=="Tiandi Grassland"] <- NA

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
NPP_all$alpha_a <- (NPP_all$alpha)
NPP_all$Tg_a <- NPP_all$Tg
NPP_all$PPFD_a <- log(NPP_all$PPFD)

NPP_all$PPFD_total_a <- log(NPP_all$PPFD_total)

NPP_all$vpd_a <- log(NPP_all$vpd)
NPP_all$fAPAR_a <- NPP_all$fAPAR
NPP_all$CNrt_a <- log(NPP_all$CNrt)
NPP_all$LMA_a <- log(NPP_all$LMA)
NPP_all$vcmax25_a <- log(NPP_all$vcmax25)

NPP_all$alpha_b <- (NPP_all$alpha_sites)
NPP_all$Tg_b <- NPP_all$Tg_sites
NPP_all$PPFD_b <- log(NPP_all$PPFD_sites)
NPP_all$vpd_b <- log(NPP_all$vpd_sites)

NPP_forest <- subset(NPP_all,pft=="Forest")

#BP_dataset <- na.omit(NPP_forest[,c("tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","alpha_a","site_a")])
#model1 <- stepwise(BP_dataset,"tnpp_a")
#model1[[1]]
#model1[[2]]
#bp_model <- (lmer(tnpp_a~Tg_a+observedfAPAR_a+obs_age_a+PPFD_a+alpha_a+(1|site_a),data=BP_dataset))
#summary(bp_model)

BP_dataset2 <- na.omit(NPP_forest[,c("tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
a2 <- stepwise(BP_dataset2,"tnpp_a")
a2[[1]]
a2[[3]]
bp_model <- (lmer(tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+(1|site_a),data=BP_dataset2))
summary(bp_model)
r.squaredGLMM(bp_model)

#try an alternative, with additional PPFD_total
#but no improvement
BP_dataset2a <- na.omit(NPP_forest[,c("tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_total_a","vpd_a","site_a")])
a2a <- stepwise(BP_dataset2a,"tnpp_a")
a2a[[1]]
a2a[[3]]
bp_modela <- (lmer(tnpp_a~Tg_a+fAPAR_a+vpd_a+CNrt_a+age_a+(1|site_a),data=BP_dataset2a))
summary(bp_modela)
r.squaredGLMM(bp_modela)

#anpp_tnpp_dataset <- na.omit(NPP_forest[,c("anpp_tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","alpha_a","site_a")])
#dim(anpp_tnpp_dataset)
#model2 <- stepwise(anpp_tnpp_dataset,"anpp_tnpp_a")
#model2[[1]]
#model2[[2]]
#anpp_tnpp_model <- (lmer(anpp_tnpp_a~soilCN_a+obs_age_a+observedfAPAR_a+(1|site_a),data=anpp_tnpp_dataset))
#summary(anpp_tnpp_model)

#mapped
anpp_tnpp_dataset2 <- na.omit(NPP_forest[,c("anpp_tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(anpp_tnpp_dataset2)
model2a <- stepwise(anpp_tnpp_dataset2,"anpp_tnpp_a")
model2a[[1]]
model2a[[3]]
anpp_tnpp_model <- (lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+age_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(anpp_tnpp_model)
r.squaredGLMM(anpp_tnpp_model)

#mapped, try an alternative, with additional PPFD_total - still negative 
#anpp_tnpp_dataset2a <- na.omit(NPP_forest[,c("anpp_tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_total_a","vpd_a","site_a")])
#dim(anpp_tnpp_dataset2a)
#model2a <- stepwise(anpp_tnpp_dataset2a,"anpp_tnpp_a")
#model2a[[1]]
#model2a[[3]]
#anpp_tnpp_model <- (lmer(anpp_tnpp_a~CNrt_a+PPFD_total_a+Tg_a+vpd_a+age_a+(1|site_a),data=anpp_tnpp_dataset2a))
#summary(anpp_tnpp_model)
#r.squaredGLMM(anpp_tnpp_model)

#only one predictor found when including measured predictors
#anpp_leafnpp_dataset <- na.omit(NPP_forest[,c("anpp_leafnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","alpha_a","site_a")])
#dim(anpp_tnpp_dataset)
#model3 <- stepwise(anpp_leafnpp_dataset,"anpp_leafnpp_a")
#model3[[1]]
#model3[[2]]
#anpp_leafnpp_model <- (lmer(anpp_leafnpp_a~obs_age_a+(1|site_a),data=anpp_leafnpp_dataset))
#summary(anpp_leafnpp_model)

# two versions - later one selected
#anpp_leafnpp_dataset <- subset(na.omit(NPP_forest[,c("anpp_leafnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","alpha_a","site_a")]))
#dim(anpp_leafnpp_dataset)
#model3 <- stepwise(anpp_leafnpp_dataset,"anpp_leafnpp_a")
#model3[[1]]
#model3[[3]]
#anpp_leafnpp_model2 <- (lmer(anpp_leafnpp_a~age_a+PPFD_a+vpd_a+(1|site_a),data=anpp_leafnpp_dataset)) 
#r.squaredGLMM(anpp_leafnpp_model2)

anpp_leafnpp_dataset <- na.omit(NPP_forest[,c("anpp_leafnpp_a","Tg_a","PPFD_a","vpd_a","site_a")])
model3 <- stepwise(anpp_leafnpp_dataset,"anpp_leafnpp_a")
model3[[1]]
model3[[3]]
anpp_leafnpp_model <- (lmer(anpp_leafnpp_a~Tg_a+vpd_a+PPFD_a+(1|site_a),data=anpp_leafnpp_dataset)) 
r.squaredGLMM(anpp_leafnpp_model)

#check tnpp grassland
#not filtering any management/non-management! while previous do so
#removing tiandi's grassland
NPP_grassland <- subset(NPP_all,pft=="Grassland" & is.na(Nmin)==TRUE)
grassland_sitemean <- aggregate(NPP_grassland,by=list(NPP_grassland$site), FUN=mean, na.rm=TRUE) 

BP_dataset_grass <- na.omit(grassland_sitemean[,c("tnpp_a","Tg_a","PPFD_a","vpd_a","CNrt_a","fAPAR_a")])
model_g1 <- stepwise_lm(BP_dataset_grass,"tnpp_a")
model_g1[[1]]
model_g1[[2]]

bp_grass_model <- (lm(tnpp_a~PPFD_a+Tg_a,data=BP_dataset_grass))
summary(bp_grass_model)
r.squaredGLMM(bp_grass_model)

#try alternative for grassland bp with total PPFD
#BP_dataset_grass2 <- na.omit(grassland_sitemean[,c("tnpp_a","Tg_a","PPFD_total_a","vpd_a","CNrt_a","fAPAR_a")])
#model_g1a <- stepwise_lm(BP_dataset_grass2,"tnpp_a")
#model_g1a[[1]]
#model_g1a[[2]]

#bp_grass_modela <- (lm(tnpp_a~vpd_a+Tg_a,data=BP_dataset_grass2))
#summary(bp_grass_modela)
#r.squaredGLMM(bp_grass_model)


#anpp/tnpp
anpp_tnpp_dataset_grass <- na.omit(grassland_sitemean[,c("anpp_tnpp_a","Tg_a","PPFD_a","vpd_a","CNrt_a","fAPAR_a")])
model_g2 <- stepwise_lm(anpp_tnpp_dataset_grass,"anpp_tnpp_a")
model_g2[[1]]
summary(lm(anpp_tnpp_a~Tg_a,data=anpp_tnpp_dataset_grass))
#non-significant! so alternatively using constant ratio
summary((lm(anpp_a~-1+tnpp_a,data=grassland_sitemean))) # 0.49 for anpp, so 0.51 for bnpp

#leaf Nmass
###2. leaf Nmass basing on a site-species model
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
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

#validation directly
sitemean$pred_nmass <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax.25/sitemean$lma
sitemean$obs_nmass <- sitemean$narea/sitemean$lma
p11 <- analyse_modobs2(sitemean,"pred_nmass","obs_nmass", type = "points")


###3. NRE model basing site-mean (lm)
NRE_climate <- read.csv("/Users/yunpeng/data/NRE_various/NRE_dataset.csv")
NRE_climate$nre_a <- log(NRE_climate$nre/(1-NRE_climate$nre))
NRE_climate$Tg_a <- NRE_climate$Tg
NRE_climate$vpd_a <- log(NRE_climate$vpd)
NRE_climate$PPFD_a <- log(NRE_climate$PPFD)
NRE_climate$alpha_a <- NRE_climate$alpha

nre_model <- lm(nre_a~Tg_a+vpd_a,data=NRE_climate)
summary(nre_model)

#validation directly
NRE_climate$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] + summary(nre_model)$coefficients[2,1] *NRE_climate$Tg_a + 
                                      summary(nre_model)$coefficients[3,1] * NRE_climate$vpd_a))))

p12 <- analyse_modobs2(NRE_climate,"pred_nre","nre", type = "points")


# Cmass constant = 46%
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
SP_input_mean <- aggregate(SP_input,by=list(SP_input$lon,SP_input$lat), FUN=mean, na.rm=TRUE)
summary(SP_input_mean$C_percent,na.rm=TRUE)
dim(subset(SP_input_mean,C_percent>0))

# Root C/N = 94 - as derived from median values of collected samples
NPP_forest_sitemean <- aggregate(NPP_forest,by=list(NPP_forest$site), FUN=mean, na.rm=TRUE) 
summary(NPP_forest_sitemean$CN_root_final) # using median = 94
dim(subset(NPP_forest,CN_root_final>0))

# Wood C/N = 100 - as derived median values of TRY database
wood_cn <- read.csv("~/data/CN_wood/wood_cn.csv")
wood_cn_sitemean <- aggregate(wood_cn,by=list(wood_cn$lon,wood_cn$lat), FUN=mean, na.rm=TRUE) 
summary(wood_cn_sitemean$OrigValueStr)
dim(wood_cn_sitemean)

##2. For Grassland

#leaf c/n model. median = 18. 
grassland_sitemean <- aggregate(NPP_grassland,by=list(NPP_grassland$site), FUN=mean, na.rm=TRUE) 
summary(grassland_sitemean$CN_leaf_final)
dim(subset(grassland_sitemean,CN_leaf_final>0))

#root c/n model median = 41.
summary(grassland_sitemean$CN_root_final)
dim(subset(grassland_sitemean,CN_root_final>0))


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
length(NRE_df$NRE)
#median of NRE in grassland is 0.69, site  =26

#check largest points
qq <- subset(NPP_grassland,TNPP_1>2000)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(qq$lon,qq$lat, col="red", pch=16,cex=1)

#final model look
#bp_model,anpp_tnpp_model,anpp_leafnpp_model,bp_grass_model,0.49/0.51,n1,nre_model

#forest validation
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
                  summary(anpp_leafnpp_model)$coefficients[2,1]* NPP_forest$Tg_a +
                  summary(anpp_leafnpp_model)$coefficients[3,1] * NPP_forest$vpd_a + 
                  summary(anpp_leafnpp_model)$coefficients[4,1] * NPP_forest$PPFD_a))))

NPP_forest$pred_wnpp <- NPP_forest$pred_anpp - NPP_forest$pred_lnpp

NPP_forest$pred_leafnc <- (summary(n1)$coefficients[1,1]/0.47) + 
  (summary(n1)$coefficients[2,1]/0.47) * NPP_forest$vcmax25/NPP_forest$LMA

NPP_forest$pred_lnf <- NPP_forest$pred_lnpp*NPP_forest$pred_leafnc
NPP_forest$pred_nre <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                                      summary(nre_model)$coefficients[2,1] *NPP_forest$Tg_a +
                                      summary(nre_model)$coefficients[3,1] * NPP_forest$vpd_a))))

NPP_forest$pred_wnf <- NPP_forest$pred_wnpp/100
NPP_forest$pred_bnf <- NPP_forest$pred_bnpp/94

NPP_forest$pred_nuptake <- NPP_forest$pred_lnf*(1-NPP_forest$pred_nre)+NPP_forest$pred_wnf+NPP_forest$pred_bnf

NPP_forest$lnf_obs_final <-NPP_forest$NPP.foliage/NPP_forest$CN_leaf_final
NPP_forest$bnf_obs_final  <- NPP_forest$BNPP_1/NPP_forest$CN_root_final
NPP_forest$wnf_obs_final  <- NPP_forest$NPP.wood/NPP_forest$CN_wood_final

#aggregate
p1 <- analyse_modobs2(NPP_forest, "pred_npp","TNPP_1",type = "points")
p2 <- analyse_modobs2(NPP_forest,"pred_anpp", "ANPP_2",type = "points")
p3 <- analyse_modobs2(NPP_forest,"pred_lnpp","NPP.foliage", type = "points")
p4 <- analyse_modobs2(NPP_forest, "pred_wnpp","NPP.wood",type = "points")
p5 <- analyse_modobs2(NPP_forest,"pred_bnpp","BNPP_1", type = "points")
p6 <- analyse_modobs2(NPP_forest,"pred_lnf","lnf_obs_final", type = "points")
p7 <- analyse_modobs2(NPP_forest,"pred_nuptake","Nmin", type = "points")

#grass validation
NPP_grassland$grassland_pred_npp <- summary(bp_grass_model)$coefficients[1,1] +  
  summary(bp_grass_model)$coefficients[2,1] * NPP_grassland$PPFD_a +
  summary(bp_grass_model)$coefficients[3,1] * NPP_grassland$Tg_a 

NPP_grassland$grassland_pred_anpp <- NPP_grassland$grassland_pred_npp * 0.49
NPP_grassland$grassland_pred_lnf <- NPP_grassland$grassland_pred_anpp/18
NPP_grassland$grassland_pred_bnpp <- NPP_grassland$grassland_pred_npp - NPP_grassland$grassland_pred_anpp
NPP_grassland$grassland_pred_bnf <- NPP_grassland$grassland_pred_bnpp/41

p8 <- analyse_modobs2(NPP_grassland,"grassland_pred_npp","TNPP_1", type = "points")
p9 <- analyse_modobs2(NPP_grassland,"grassland_pred_anpp","ANPP_2", type = "points")
p10 <- analyse_modobs2(NPP_grassland,"grassland_pred_bnpp","BNPP_1", type = "points")
white <- theme(plot.background=element_rect(fill="white", color="white"))

#fig.2 validation
plot_grid(p1$gg,p2$gg,p5$gg,
          p8$gg,p9$gg,p10$gg,
          p3$gg,p4$gg,p11$gg,
          p6$gg,p12$gg,p7$gg, 
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'),
          ncol=3,label_x = 0.9,label_y=0.92)+white
ggsave(paste("~/data/output/newphy_fig2.jpg",sep=""),width = 20, height = 20)

#now, inputting all predictors
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

alpha <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/alpha.nc"),
  varnam = "alpha"))

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

###input land cover
ncin <- nc_open("~/data/landcover/modis_landcover_halfdeg_2010_FILLED.nc")
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon) 
lat<-ncvar_get(ncin,"lat")
nlat<-dim(lat)
pftcover <-ncvar_get(ncin,"pftcover")
nc_close(ncin)
pftcover_long <- as.vector(pftcover)
pftcover <- as.data.frame(matrix(pftcover_long, nrow = nlon * nlat, ncol = 10))
#see get_fpc_grid function: https://github.com/stineb/sofun/blob/db7a9e8e486f576fd7b9f1f74edb1df7a8d2c4f7/src/forcing_global_wmodel.mod.f90 
#it clarified that: 1-6 is forest, 8 is grassland
forest_percent <- rowSums(pftcover[,1:6],na.rm=TRUE)/rowSums(pftcover[,c(1:6,8)],na.rm = TRUE)
grass_percent <- pftcover[,8]/rowSums(pftcover[,c(1:6,8)],na.rm = TRUE)
summary(grass_percent + forest_percent) # check - their sum = 1, perfect!

###3. calculate weighted-sum
#firstly - filter na points - so that all output map has same numbers of NA.
all_predictors <- as.data.frame(cbind(Tg$myvar,PPFD$myvar,vpd$myvar,
                                      fAPAR$myvar,age$myvar,alpha$myvar,
                                      CNrt$myvar,LMA$myvar,vcmax25_df$myvar))
all_predictors$available_grid = rowMeans(all_predictors)
#just to find all na columns
all_predictors$available_grid[is.na(all_predictors$available_grid)==FALSE] <- 1
summary(all_predictors$available_grid)
available_grid2 <- all_predictors$available_grid

#represent grids when stand-age is especially in NA, but others are fine
names(all_predictors) <- c("Tg","PPFD","vpd","fAPAR","age","alpha","CNrt","LMA","vcmax25","available_grid")
all_predictors$lon <- vcmax25_df$lon
all_predictors$lat <- vcmax25_df$lat
summary(all_predictors)

#final calculation - now divide into forest, grassland and pft
#available_grid2 here was used as a list of data to identify if a grid is available (=1) or any prediction fields shown as NA 

#bp_model,anpp_tnpp_model,anpp_leafnpp_model,bp_grass_model,0.49/0.51,n1,nre_model

npp_f <- summary(bp_model)$coefficients[1,1] +  
  summary(bp_model)$coefficients[2,1] * Tg$myvar +
  summary(bp_model)$coefficients[3,1] * fAPAR$myvar +
  summary(bp_model)$coefficients[4,1] * log(PPFD$myvar) +
  summary(bp_model)$coefficients[5,1] * log(CNrt$myvar)+
  summary(bp_model)$coefficients[6,1] * log(age$myvar)

#show npp<=0 points - high lat!
#newmap <- getMap(resolution = "low")
#plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
#points(subset(as.data.frame(cbind(vcmax25_df,npp_f)),npp_f<=0)$lon,subset(as.data.frame(cbind(vcmax25_df,npp_f)),npp_f<=0)$lat, col="red", pch=16,cex=1)

#convert to NA!
npp_f[npp_f<=0] <- NA

anpp_f <- npp_f * 
  (1/(1 + exp(-(summary(anpp_tnpp_model)$coefficients[1,1]+
                  summary(anpp_tnpp_model)$coefficients[2,1] * log(CNrt$myvar) +
                  summary(anpp_tnpp_model)$coefficients[3,1] * log(PPFD$myvar) + 
                  summary(anpp_tnpp_model)$coefficients[4,1] * Tg$myvar+
                  summary(anpp_tnpp_model)$coefficients[5,1] * log(age$myvar)))))

bnpp_f <- npp_f - anpp_f

lnpp_f <- anpp_f * (1/(1 + exp(-(summary(anpp_leafnpp_model)$coefficients[1,1]+
                  summary(anpp_leafnpp_model)$coefficients[2,1]* Tg$myvar +
                  summary(anpp_leafnpp_model)$coefficients[3,1] * vpd$myvar + 
                  summary(anpp_leafnpp_model)$coefficients[4,1] * log(PPFD$myvar)))))

wnpp_f <- anpp_f - lnpp_f

leafnc_f <- (summary(n1)$coefficients[1,1]/0.47) + 
  (summary(n1)$coefficients[2,1]/0.47) * vcmax25_df$myvar/LMA$myvar

nre_f <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                       summary(nre_model)$coefficients[2,1] *Tg$myvar +
                       summary(nre_model)$coefficients[3,1] * log(vpd$myvar)))))

lnf_f <- (1-nre_f)* leafnc_f * lnpp_f
wnf_f <- wnpp_f/100
#100 is constant wood c/n
bnf_f <- bnpp_f/94
#94 is constant root c/n
nuptake_f <- lnf_f + wnf_f + bnf_f

#grass
npp_g <- summary(bp_grass_model)$coefficients[1,1] +  
  summary(bp_grass_model)$coefficients[2,1] * log(PPFD$myvar) +
  summary(bp_grass_model)$coefficients[3,1] * Tg$myvar 

#show npp<=0 points - high lat!
#newmap <- getMap(resolution = "low")
#plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
#points(subset(as.data.frame(cbind(vcmax25_df,npp_g)),npp_g<=0)$lon,subset(as.data.frame(cbind(vcmax25_df,npp_g)),npp_g<=0)$lat, col="red", pch=16,cex=1)

#convert to NA!
npp_g[npp_g<=0] <- NA

summary(npp_g)

anpp_g <- npp_g* 0.49
bnpp_g <- npp_g-anpp_g

leafnc_g <- 1/18

nre_g <- 0.69

lnf_g <- anpp_g*leafnc_g*(1-nre_g)

bnf_g <- bnpp_g *(1/41)
#41 is constant root c/n
nuptake_g <- lnf_g + bnf_g

#combine into 2 pfts
npp_pft <- available_grid2* (npp_f*forest_percent +npp_g*grass_percent)
npp_forest <- available_grid2* (npp_f*forest_percent)
npp_grass <- available_grid2* (npp_g*grass_percent)

anpp_pft <- available_grid2*(anpp_f*forest_percent +anpp_g*grass_percent)
anpp_forest <- available_grid2* (anpp_f*forest_percent)
anpp_grass <- available_grid2* (anpp_g*grass_percent)

lnpp_forest <- available_grid2*lnpp_f*forest_percent

wnpp_forest <- available_grid2*wnpp_f*forest_percent

bnpp_pft <- available_grid2*(bnpp_f*forest_percent +bnpp_g*grass_percent)
bnpp_forest <- available_grid2* (bnpp_f*forest_percent)
bnpp_grass <- available_grid2* (bnpp_g*grass_percent)

leafcn_pft <- 1/(available_grid2*(leafnc_f*forest_percent +leafnc_g*grass_percent))
leafcn_forest <- 1/available_grid2*leafnc_f*forest_percent
leafcn_grassland <- 1/available_grid2*leafnc_g*grass_percent
summary(leafcn_pft)

nre_pft <- available_grid2*(nre_f*forest_percent +nre_g*grass_percent)
nre_forest <-  available_grid2*nre_f*forest_percent
nre_grassland <- available_grid2*nre_g*grass_percent
summary(nre_pft)

lnf_pft <- available_grid2*(lnf_f*forest_percent +lnf_g*grass_percent)
lnf_forest <- available_grid2* (lnf_f*forest_percent)
lnf_grass <- available_grid2* (lnf_g*grass_percent)

wnf_forest <- available_grid2*wnf_f*forest_percent

bnf_pft <- available_grid2*(bnf_f*forest_percent +bnf_g*grass_percent)
bnf_forest <- available_grid2* (bnf_f*forest_percent)
bnf_grass <- available_grid2* (bnf_g*grass_percent)

nuptake_pft <- available_grid2*(nuptake_f*forest_percent +nuptake_g*grass_percent)
nuptake_pft_final <- nuptake_pft
nuptake_forest <- available_grid2* (nuptake_f*forest_percent)
nuptake_grass <- available_grid2* (nuptake_g*grass_percent)

all_maps <- as.data.frame(cbind(vcmax25_df,npp_pft,npp_forest,npp_grass,
                                anpp_pft,anpp_forest,anpp_grass,
                                bnpp_pft,bnpp_forest,bnpp_grass,
                                lnpp_forest,wnpp_forest,wnf_forest,
                                leafcn_pft,leafcn_forest,leafcn_grassland,
                                nre_pft,nre_forest,nre_grassland,
                                lnf_pft,lnf_forest,lnf_grass,
                                bnf_pft,bnf_forest,bnf_grass,
                                nuptake_pft,nuptake_forest,nuptake_grass))

summary(all_maps)

#####area_m2 to show each grid's area in m2
calc_area <- function( lat, dx=1, dy=1 ){
  r_earth <- 6370499.317638  # to be consistent with how Ferret calculates areas of spheres (https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2016/msg00155.html)
  area <- 4 * r_earth^2 * 0.5 * dx * pi/180 * cos( abs(lat) * pi/180 ) * sin( 0.5 * dy * pi/180 )
  return(area)
}
lonlat <- vcmax25_df[,c("lon","lat")]
area_m2 <- calc_area(lonlat$lat,0.5,0.5)
#fland - to show each grid's land cover percentage
nc <- read_nc_onefile("~/data/fland/global.fland.nc") #Input nc
output_fland <- nc_to_df(nc, varnam = "fland")
fland <- output_fland$myvar
#include conversion factor (from g to Pg)
conversion <- area_m2 * fland /1e+15

#Fig.3 global maps
white <- theme(plot.background=element_rect(fill="white", color="white"))
npp_f1 <- (NPP_forest %>% filter(TNPP_1>0) %>% filter(pred_npp>0))[,c("lon","lat")]
npp_g1 <- (NPP_grassland %>% filter(TNPP_1>0) %>% filter(grassland_pred_npp>0))[,c("lon","lat")]
total_value <- round(sum(all_maps[,"npp_pft"]*conversion,na.rm=TRUE),2)

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","npp_pft")]),
                varnam = "npp_pft",latmin = -65, latmax = 85,combine=FALSE)

a3 <- gg$ggmap +
  geom_point(data=npp_f1,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=npp_g1,aes(lon,lat),col="blue",size=1.5)+
  labs(title = paste("BP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white

a4 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))+ white

#3.anpp
anpp_f1 <- (NPP_forest %>% filter(ANPP_2>0) %>% filter(pred_anpp>0))[,c("lon","lat")]
anpp_g1 <- (NPP_grassland %>% filter(ANPP_2>0) %>% filter(grassland_pred_anpp>0))[,c("lon","lat")]
total_value <- round(sum(all_maps[,"anpp_pft"]*conversion,na.rm=TRUE),2)

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","anpp_pft")]),
                varnam = "anpp_pft",latmin = -65, latmax = 85,combine=FALSE)

a5 <- gg$ggmap +
  geom_point(data=anpp_f1,aes(lon,lat),col="red",size=1.5)+
  geom_point(data=anpp_g1,aes(lon,lat),col="blue",size=1.5)+
  labs(title = paste("ANPP:", total_value, "PgC/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white

a6 <- gg$gglegend+labs(title = ~paste("gC m"^-2,"yr"^-1))+ white

#4. leaf c/n
SP_input <- read.csv(file="~/data/leaf_traits/combined_leaf_traits_updated.csv") #new one 
SP_input2 <- SP_input[,c("lat","lon","z","Vcmax25","narea","lma")]
sitemean <- aggregate(SP_input2,by=list(SP_input2$lon,SP_input2$lat), FUN=mean, na.rm=TRUE) 
sitemean$pred_leafn <- (summary(n1)$coefficients[1,1]) + (summary(n1)$coefficients[2,1])* sitemean$Vcmax25/sitemean$lma
sitemean$obs_leafn <- sitemean$narea/sitemean$lma

laefcn <- (sitemean %>% filter(pred_leafn>0) %>% filter(obs_leafn>0))[,c("lon","lat")]

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","leafcn_pft")]),
                varnam = "leafcn_pft",latmin = -65, latmax = 85,combine=FALSE)

total_value <- round(mean(leafcn_pft,na.rm=TRUE),2)
a7 <- gg$ggmap +
  geom_point(data=laefcn,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("Leaf C/N: ", total_value))+
  theme_grey(base_size = 12)+ white

a8 <- gg$gglegend+ white

#5. NRE
NRE_validation <- read.csv(file="~/data/NPP_final/NRE_validation.csv") 
nre_site <- (NRE_validation %>% filter(pred_nre>0) %>% filter(NRE>0))[,c("lon","lat")]

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nre_pft")]),
                varnam = "nre_pft",latmin = -65, latmax = 85,combine=FALSE)

total_value <- round(mean(nre_pft,na.rm=TRUE),2)

a9 <- gg$ggmap +
  geom_point(data=nre_site,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("NRE: ", total_value))+
  theme_grey(base_size = 12)+ white

a10 <- gg$gglegend+ white

#6. nuptake
nuptake_f1 <- (NPP_forest %>% filter(pred_nuptake>0) %>% filter(Nmin>0))[,c("lon","lat")]

total_value <- 1000*round(sum(all_maps[,"nuptake_pft"]*conversion,na.rm=TRUE),2) #unit convert from PgN/yr to TgN/yr

gg <- plot_map3(na.omit(all_maps[,c("lon","lat","nuptake_pft")]),
                varnam = "nuptake_pft",latmin = -65, latmax = 85,combine=FALSE)

a11 <- gg$ggmap +
  geom_point(data=nuptake_f1,aes(lon,lat),col="red",size=1.5)+
  labs(title = paste("N uptake: ", total_value, "TgN/yr", sep=" " ))+
  theme_grey(base_size = 12)+ white

a12 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))+ white

#NUE
all_maps$NUE <- all_maps$npp_pft/all_maps$nuptake_pft
summary(all_maps$NUE)
gg <- plot_map3(na.omit(all_maps[,c("lon","lat","NUE")]),
                varnam = "NUE",latmin = -65, latmax = 85,combine=FALSE)

total_value<- round(sum(all_maps[,"npp_pft"]*conversion,na.rm=TRUE)/sum(all_maps[,"nuptake_pft"]*conversion,na.rm=TRUE),2)

a13 <- gg$ggmap +
  labs(title = paste("NUE: ", total_value, " ", sep=" " ))+
  theme_grey(base_size = 12)+ white

a14 <- gg$gglegend+labs(title = ~paste("gN m"^-2,"yr"^-1))+ white

plot_grid(a3,a4,a5,a6,a7,a8,
          a9,a10,a11,a12,a13,a14,
          nrow=2,
          rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' '))

ggsave(paste("~/data/output/newphy_fig3.jpg",sep=""),width = 20, height = 10*(2/3))

#work on effect of each factor on NUE
#now, work on NUE
nue_final <- all_maps$NUE # by default
cal_nue <- function(Tg_pred,PPFD_pred,vpd_pred,fAPAR_pred,age_pred,CNrt_pred,LMA_pred,vcmax25_pred){
  npp_f <- summary(bp_model)$coefficients[1,1] +  
    summary(bp_model)$coefficients[2,1] * Tg_pred +
    summary(bp_model)$coefficients[3,1] * fAPAR_pred +
    summary(bp_model)$coefficients[4,1] * log(PPFD_pred) +
    summary(bp_model)$coefficients[5,1] * log(CNrt_pred)+
    summary(bp_model)$coefficients[6,1] * log(age_pred)
  
  npp_f[npp_f<=0] <- NA
  
  anpp_f <- npp_f * 
    (1/(1 + exp(-(summary(anpp_tnpp_model)$coefficients[1,1]+
                    summary(anpp_tnpp_model)$coefficients[2,1] * log(CNrt_pred) +
                    summary(anpp_tnpp_model)$coefficients[3,1] * log(PPFD_pred) + 
                    summary(anpp_tnpp_model)$coefficients[4,1] * Tg_pred+
                    summary(anpp_tnpp_model)$coefficients[5,1] * log(age_pred)))))
  
  bnpp_f <- npp_f - anpp_f
  
  lnpp_f <- anpp_f * (1/(1 + exp(-(summary(anpp_leafnpp_model)$coefficients[1,1]+
                                     summary(anpp_leafnpp_model)$coefficients[2,1]* Tg_pred +
                                     summary(anpp_leafnpp_model)$coefficients[3,1] * vpd_pred + 
                                     summary(anpp_leafnpp_model)$coefficients[4,1] * log(PPFD_pred)))))
  
  wnpp_f <- anpp_f - lnpp_f
  
  leafnc_f <- (summary(n1)$coefficients[1,1]/0.47) + 
    (summary(n1)$coefficients[2,1]/0.47) * vcmax25_pred/LMA_pred
  
  nre_f <- (1/(1+exp(-(summary(nre_model)$coefficients[1,1] +
                         summary(nre_model)$coefficients[2,1] *Tg_pred +
                         summary(nre_model)$coefficients[3,1] * log(vpd_pred)))))
  
  lnf_f <- (1-nre_f)* leafnc_f * lnpp_f
  wnf_f <- wnpp_f/100
  #100 is constant wood c/n
  bnf_f <- bnpp_f/94
  #94 is constant root c/n
  nuptake_f <- lnf_f + wnf_f + bnf_f
  
  #grass
  npp_g <- summary(bp_grass_model)$coefficients[1,1] +  
    summary(bp_grass_model)$coefficients[2,1] * log(PPFD_pred) +
    summary(bp_grass_model)$coefficients[3,1] * Tg_pred 
  
  anpp_g <- npp_g* 0.49
  bnpp_g <- npp_g-anpp_g
  
  leafnc_g <- 1/18
  
  nre_g <- 0.69
  
  lnf_g <- anpp_g*leafnc_g*(1-nre_g)
  
  bnf_g <- bnpp_g *(1/41)
  #41 is constant root c/n
  nuptake_g <- lnf_g + bnf_g
  
  #combine into 2 pfts
  npp_pft <- available_grid2* (npp_f*forest_percent +npp_g*grass_percent)
  npp_forest <- available_grid2* (npp_f*forest_percent)
  npp_grass <- available_grid2* (npp_g*grass_percent)
  
  anpp_pft <- available_grid2*(anpp_f*forest_percent +anpp_g*grass_percent)
  anpp_forest <- available_grid2* (anpp_f*forest_percent)
  anpp_grass <- available_grid2* (anpp_g*grass_percent)
  
  lnpp_forest <- available_grid2*lnpp_f*forest_percent
  
  wnpp_forest <- available_grid2*wnpp_f*forest_percent
  
  bnpp_pft <- available_grid2*(bnpp_f*forest_percent +bnpp_g*grass_percent)
  bnpp_forest <- available_grid2* (bnpp_f*forest_percent)
  bnpp_grass <- available_grid2* (bnpp_g*grass_percent)
  
  leafcn_pft <- 1/(available_grid2*(leafnc_f*forest_percent +leafnc_g*grass_percent))
  leafcn_forest <- 1/available_grid2*leafnc_f*forest_percent
  leafcn_grassland <- 1/available_grid2*leafnc_g*grass_percent
  summary(leafcn_pft)
  
  nre_pft <- available_grid2*(nre_f*forest_percent +nre_g*grass_percent)
  nre_forest <-  available_grid2*nre_f*forest_percent
  nre_grassland <- available_grid2*nre_g*grass_percent
  summary(nre_pft)
  
  lnf_pft <- available_grid2*(lnf_f*forest_percent +lnf_g*grass_percent)
  lnf_forest <- available_grid2* (lnf_f*forest_percent)
  lnf_grass <- available_grid2* (lnf_g*grass_percent)
  
  wnf_forest <- available_grid2*wnf_f*forest_percent
  
  bnf_pft <- available_grid2*(bnf_f*forest_percent +bnf_g*grass_percent)
  bnf_forest <- available_grid2* (bnf_f*forest_percent)
  bnf_grass <- available_grid2* (bnf_g*grass_percent)
  
  nuptake_pft <- available_grid2*(nuptake_f*forest_percent +nuptake_g*grass_percent)
  nuptake_forest <- available_grid2* (nuptake_f*forest_percent)
  nuptake_grass <- available_grid2* (nuptake_g*grass_percent)
  
  new_nue <-npp_pft/nuptake_pft
  
  nue_pft_ratio <- log(nue_final/new_nue)
  
  return(nue_pft_ratio)
}

nue_standard <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,
                        age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
summary(nue_standard)

mean_Tg <- mean(Tg$myvar[is.na(available_grid2)==FALSE]);mean_Tg
mean_PPFD <- mean(PPFD$myvar[is.na(available_grid2)==FALSE]);mean_PPFD
mean_vpd <- mean(vpd$myvar[is.na(available_grid2)==FALSE]);mean_vpd
mean_fAPAR <- mean(fAPAR$myvar[is.na(available_grid2)==FALSE]);mean_fAPAR
mean_age <- mean(age$myvar[is.na(available_grid2)==FALSE]);mean_age
mean_CNrt <- mean(CNrt$myvar[is.na(available_grid2)==FALSE]);mean_CNrt
mean_LMA <- mean(LMA$myvar[is.na(available_grid2)==FALSE]);mean_LMA
mean_vcmax25 <- mean(vcmax25_df$myvar[is.na(available_grid2)==FALSE]);mean_vcmax25
mean_alpha <- mean(alpha$myvar[is.na(available_grid2)==FALSE]);mean_alpha


nue_Tg <- cal_nue(rep(mean_Tg,259200),PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_PPFD <- cal_nue(Tg$myvar,rep(mean_PPFD,259200),vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_vpd <- cal_nue(Tg$myvar,PPFD$myvar,rep(mean_vpd,259200),fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_fAPAR <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,rep(mean_fAPAR,259200),age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_age <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,rep(mean_age,259200),CNrt$myvar,LMA$myvar,vcmax25_df$myvar)
nue_CNrt <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,rep(mean_CNrt,259200),LMA$myvar,vcmax25_df$myvar)
nue_LMA <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,rep(mean_LMA,259200),vcmax25_df$myvar)
nue_vcmax25 <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,rep(mean_vcmax25,259200))
#nue_alpha <- cal_nue(Tg$myvar,PPFD$myvar,vpd$myvar,fAPAR$myvar,age$myvar,CNrt$myvar,LMA$myvar,vcmax25_df$myvar)

nue_all <- as.data.frame(cbind(all_predictors,nue_Tg,nue_PPFD,nue_vpd,nue_fAPAR,nue_age,nue_CNrt,nue_LMA,nue_vcmax25))
summary(nue_all)
nue_all$NUE <- all_maps$NUE
#show soil C/N effects on NUE
#aggregate
nue_all$name[nue_all$lat< -60] <- "90°S ~ 60°S"
nue_all$name[nue_all$lat >= -60 & nue_all$lat < -30] <- "60°S ~ 30°S"
nue_all$name[nue_all$lat >= -30 & nue_all$lat < 0] <- "30°S ~ 0"
nue_all$name[nue_all$lat >= 0 & nue_all$lat < 30]<- "0 ~ 30°N"
nue_all$name[nue_all$lat >= 30 & nue_all$lat < 60]<- "30°N ~ 60°N"
nue_all$name[nue_all$lat >= 60]<- "60°N ~ 90°N"

nue_all$name <- factor(nue_all$name,levels = c("90°S ~ 60°S","60°S ~ 30°S","30°S ~ 0",
                                               "0 ~ 30°N","30°N ~ 60°N","60°N ~ 90°N"))
nue_all <- subset(nue_all,name!="90°S ~ 60°S")
a1 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_CNrt),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Soil C/N effect") 
a2 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_age),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Age effect") 
a3 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_fAPAR),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "fAPAR effect") 
a4 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_Tg),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Tg effect") 
a5 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_vpd),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "vpd effect") 
a6 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_PPFD),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "PPFD effect") 
a7 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_vcmax25),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "Vcmax25 effect") 
a8 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = nue_LMA),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")+ geom_vline(xintercept=0, linetype="dashed")+labs(y= " ", x = "LMA effect") 
a9 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = NUE),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="blue")

aa1 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = CNrt),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Soil C/N") 
aa2 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = age),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Age") 
aa3 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = fAPAR),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "fAPAR") 
aa4 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = Tg),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Tg") 
aa5 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = vpd),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "vpd") 
aa6 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = PPFD),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "PPFD") 
aa7 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = vcmax25),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "Vcmax25") 
aa8 <- ggplot(data = nue_all) + stat_summary(mapping = aes(y = name, x = LMA),fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median,color="black")+labs(y= " ", x = "LMA") 

plot_grid(a1,a2,a3,a4,a5,a6,a7,a8,a9)+white
ggsave(paste("~/data/output/newphy_fig4.jpg",sep=""),width = 20, height = 10)

plot_grid(aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8)+white
ggsave(paste("~/data/output/newphy_fig4_parallel.jpg",sep=""),width = 20, height = 10)

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_CNrt")]),varnam = "nue_CNrt",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
g1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15);g2 <- gg$gglegend

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_age")]),varnam = "nue_age",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4","royalblue3","wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
g3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15);g4 <- gg$gglegend+labs(title = ~paste("years"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_fAPAR")]),varnam = "nue_fAPAR",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4","royalblue3", "wheat","tomato3"),
                breaks = seq(-0.04,0.04,0.01))
g5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15);g6 <- gg$gglegend

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_PPFD")]),varnam = "nue_PPFD",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3","tomato4"),
                breaks = seq(-0.08,0.12,0.02))
g7 <- gg$ggmap +labs(title = "PPFD")+theme_grey(base_size = 15);g8 <- gg$gglegend + labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_Tg")]),varnam = "nue_Tg",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.30,0.30,0.10))
g9 <- gg$ggmap +labs(title =~paste(T[g]))+theme_grey(base_size = 15);g10 <- gg$gglegend + labs(title =~paste("\u00B0C"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_vpd")]),varnam = "nue_vpd",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.20,0.30,0.05))
g11 <- gg$ggmap +labs(title = "D")+theme_grey(base_size = 15);g12 <- gg$gglegend+ labs(title =~paste("kPa"))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_vcmax25")]), varnam = "nue_vcmax25",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.2,0.2,0.04))
g15 <- gg$ggmap +labs(title =~paste(V[cmax25]))+theme_grey(base_size = 15);g16 <- gg$gglegend+ labs(title =~paste(mu, "mol m"^-2,"s"^-1))

gg <- plot_map3(na.omit(nue_all[,c("lon","lat","nue_LMA")]),varnam = "nue_LMA",latmin = -65, latmax = 85,combine=FALSE,
                colorscale = c( "royalblue4", "wheat","tomato3"),
                breaks = seq(-0.05,0.05,0.01))
g17 <- gg$ggmap +labs(title = "LMA")+theme_grey(base_size = 15);g18 <- gg$gglegend+ labs(title =~paste("g ","m"^-2))

plot_grid(g1,g2,g3,g4,g5,g6,
          g7,g8,g9,g10,g11,g12,
          g15,g16,g17,g18,
          nrow=3,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' ',
                     '(g)',' ','(h)',' '),label_size = 15)+white

ggsave(paste("~/data/output/newphy_fig4_not_used.jpg",sep=""),width = 20, height = 10)


#figs2 representing all predictors
gg <- plot_map3(na.omit(CNrt[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d1 <- gg$ggmap +labs(title = "Soil C/N")+theme_grey(base_size = 15)+ white;d2 <- gg$gglegend+ white

gg <- plot_map3(na.omit(age[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d3 <- gg$ggmap +labs(title = "Age")+theme_grey(base_size = 15)+ white;d4 <- gg$gglegend+labs(title = ~paste("years"))+ white

gg <- plot_map3(na.omit(fAPAR[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d5 <- gg$ggmap +labs(title = "fAPAR")+theme_grey(base_size = 15)+ white;d6 <- gg$gglegend+ white

gg <- plot_map3(na.omit(PPFD[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE,breaks = seq(100,800,100))
d7 <- gg$ggmap +labs(title = "PPFD")+theme_grey(base_size = 15)+ white;d8 <- gg$gglegend + labs(title =~paste(mu, "mol m"^-2,"s"^-1))+ white

gg <- plot_map3(na.omit(Tg[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d9 <- gg$ggmap +labs(title =~paste(T[g]))+theme_grey(base_size = 15)+ white;d10 <- gg$gglegend + labs(title =~paste("\u00B0C"))+ white

gg <- plot_map3(na.omit(vpd[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d11 <- gg$ggmap +labs(title = "D")+theme_grey(base_size = 15)+ white;d12 <- gg$gglegend+ labs(title =~paste("kPa"))+ white

gg <- plot_map3(na.omit(vcmax25_df[,c("lon","lat","myvar")]), varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d13 <- gg$ggmap +labs(title =~paste(V[cmax25]))+theme_grey(base_size = 15)+ white;d14 <- gg$gglegend+ labs(title =~paste(mu, "mol m"^-2,"s"^-1))+ white

gg <- plot_map3(na.omit(LMA[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85,combine=FALSE)
d15 <- gg$ggmap +labs(title = "LMA")+theme_grey(base_size = 15)+ white;d16 <- gg$gglegend+ labs(title =~paste("g ","m"^-2))+ white

plot_grid(d1,d2,d3,d4,d5,d6,
          d7,d8,d9,d10,d11,d12,
          d13,d14,d15,d16,
          nrow=3,rel_widths = c(3/12, 1/12,3/12,1/12,3/12,1/12),
          labels = c('(a)',' ','(b)',' ','(c)',' ',
                     '(d)',' ','(e)',' ','(f)',' ',
                     '(g)',' ','(h)',' '),label_size = 15)+ white

ggsave(paste("~/data/output/newphy_figS2.jpg",sep=""),width = 20, height = 10)

#trendy
#model output
#CABLE-POP
CABLE_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CABLE-POP_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CABLE_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CABLE-POP_S2_npp_ANN_mean.nc"), varnam = "npp"))

#CABLE-POP
CLASS_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CLASS-CTEM_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CLASS_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CLASS-CTEM_S2_npp_ANN_mean.nc"), varnam = "npp"))

#CLM
#CLM_fNup <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CLM5.0_S2_fNup_ANN_mean.nc"), varnam = "fNup")) 
#this product is confusing (1) unit needs /1000 to get gn/m2/year? (2)still many values very high
CLM_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CLM5.0_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CLM_GPP$lon[CLM_GPP$lon>180] <- CLM_GPP$lon[CLM_GPP$lon>180]-360
CLM_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CLM5.0_S2_npp_ANN_mean.nc"), varnam = "npp"))
CLM_NPP$lon[CLM_NPP$lon>180] <- CLM_NPP$lon[CLM_NPP$lon>180]-360

#ISAM
#(given in kgC m-2 month-1 - they might not need to multiply with 12
#cdo seltimestep,127/156 ISAM_S2_fNup.nc a1.nc
#cdo -O timmean a1.nc ISAM_S2_fNup_ANN_mean.nc

ISAM_fNup <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ISAM_S2_fNup_ANN_mean.nc"), varnam = "fNup"))
ISAM_fNup$myvar<- ISAM_fNup$myvar*1000*31556952

ISAM_gpp <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ISAM_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
ISAM_gpp$myvar<- ISAM_gpp$myvar*1000*31556952

ISAM_npp <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ISAM_S2_npp_ANN_mean.nc"), varnam = "npp"))
ISAM_npp$myvar<- ISAM_npp$myvar*1000*31556952

#ISBA
ncin <- nc_open("/Users/yunpeng/data/trendy/v8/ISBA-CTRIP_S2_gpp_ANN_mean.nc")
lon <- ncvar_get(ncin,"lon_FULL");nlon <- dim(lon) 
lat <- ncvar_get(ncin,"lat_FULL");nlat <- dim(lat) 
ISBA_GPP <- ncvar_get(ncin,"gpp")
nc_close(ncin)
ISBA_GPP <- as.vector(ISBA_GPP)
lonlat <- expand.grid(lon,lat)
ISBA_GPP <- as.data.frame(cbind(lonlat,ISBA_GPP))
names(ISBA_GPP) <- c("lon","lat","GPP")

ncin <- nc_open("/Users/yunpeng/data/trendy/v8/ISBA-CTRIP_S2_npp_ANN_mean.nc")
ISBA_NPP <- ncvar_get(ncin,"npp")
nc_close(ncin)
ISBA_NPP <- as.vector(ISBA_NPP)
ISBA_NPP <- as.data.frame(cbind(lonlat,ISBA_NPP))
names(ISBA_NPP) <- c("lon","lat","NPP")

#JSBACH
JSBACH_fNup <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/JSBACH_S2_fNup_ANN_mean.nc"), varnam = "fNup"))
JSBACH_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/JSBACH_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
JSBACH_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/JSBACH_S2_npp_ANN_mean.nc"), varnam = "npp"))
plot_map3(na.omit(JSBACH_fNup[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85)

#JULES
# data from output/JULES-ES-1.0/JULES-ES.1p0.vn5.4.50.CRUJRA2.TRENDYv8.365.S2_Monthly_npp.nc
JULES_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/jule_gpp.nc"), varnam = "gpp"))
JULES_GPP$lon[JULES_GPP$lon>180] <- JULES_GPP$lon[JULES_GPP$lon>180]-360
JULES_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/jule_npp.nc"), varnam = "npp"))
JULES_NPP$lon[JULES_NPP$lon>180] <- JULES_NPP$lon[JULES_NPP$lon>180]-360
#looks ok 
#plot_map3(na.omit(JULES_GPP[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85)

#LPJ
LPJ_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/LPJ-GUESS_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
LPJ_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/LPJ-GUESS_S2_npp_ANN_mean.nc"), varnam = "npp"))

#ORCHIDEE
ORCHIDEE_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ORCHIDEE_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
ORCHIDEE_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ORCHIDEE_S2_npp_ANN_mean.nc"), varnam = "npp"))

#ORCHIDEE-CNP
nc_flip_lat <- function(nc){
  
  nc$lat <- rev(nc$lat)
  
  # nlat <- length(nc$lat)
  # nc$vars[[1]] <- nc$vars[[1]][,nlat:1]
  
  arr_flip_lat <- function(arr){
    nlat <- dim(arr)[2]
    arr <- arr[,nlat:1]
    return(arr)
  }
  nc$vars <- purrr::map(nc$vars[1], ~arr_flip_lat(.))
  
  return(nc)
}
ORCHICNP_fNup <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ORCHIDEE-CNP_S2_fNup_ANN_mean.nc"), varnam = "fNup"))
ORCHICNP_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ORCHIDEE-CNP_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
ORCHICNP_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ORCHIDEE-CNP_S2_npp_ANN_mean.nc"), varnam = "npp"))

#SDGVM
SDGVM_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/SDGVM_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
SDGVM_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/SDGVM_S2_npp_ANN_mean.nc"), varnam = "npp"))

#all_maps
not_used <- list(CLASS_GPP,CLASS_NPP,CLM_GPP,CLM_NPP,JSBACH_fNup,JSBACH_GPP,JSBACH_NPP)
allmaps<- list(CABLE_GPP,CABLE_NPP,ISAM_fNup,ISAM_gpp,ISAM_npp,ISBA_GPP,ISBA_NPP,
               JULES_GPP,JULES_NPP,LPJ_GPP,LPJ_NPP,ORCHIDEE_GPP,ORCHIDEE_NPP,
               ORCHICNP_fNup,ORCHICNP_GPP,ORCHICNP_NPP,SDGVM_GPP,SDGVM_NPP)
obj_name <- c("CABLE_GPP","CABLE_NPP","ISAM_fNup","ISAM_gpp","ISAM_npp","ISBA_GPP","ISBA_NPP",
              "JULES_GPP","JULES_NPP","LPJ_GPP","LPJ_NPP","ORCHIDEE_GPP","ORCHIDEE_NPP",
              "ORCHICNP_fNup","ORCHICNP_GPP","ORCHICNP_NPP","SDGVM_GPP","SDGVM_NPP")
#aggregate based on lon and lat firstly
sitemean <- unique(NPP_forest[,c("lon","lat")])
sp_sites <- SpatialPoints(sitemean) # only select lon and lat

sitemean_final <- unique(NPP_forest[,c("lon","lat")])

for (i in c(1:length(allmaps))){
  df1 <- allmaps[[i]]
  names(df1) <- c("lon","lat",obj_name[i])
  coordinates(df1) <- ~lon+lat 
  gridded(df1) <- TRUE
  df1_global <- raster(df1, obj_name[i]) 
  
  sp_sites_new <- raster::extract(df1_global, sp_sites, sp = TRUE) %>% as_tibble() %>% 
    right_join(sitemean, by = c("lon", "lat"))
  sitemean_final[,i+2] <- sp_sites_new[,1]
}

summary(sitemean_final) #convert 0 to NA
sitemean_final$CABLE_NPP[sitemean_final$CABLE_NPP==0] <- NA
sitemean_final$ISAM_fNup[sitemean_final$ISAM_fNup==0] <- NA
sitemean_final$ISAM_gpp[sitemean_final$ISAM_gpp==0] <- NA
sitemean_final$ISAM_npp[sitemean_final$ISAM_npp==0] <- NA
sitemean_final$JULES_GPP[sitemean_final$JULES_GPP==0] <- NA
sitemean_final$JULES_NPP[sitemean_final$JULES_NPP==0] <- NA

#first, our model
bp1a <- visreg(bp_model,"Tg_a",type="contrast");bp1b <- visreg(bp_model,"fAPAR_a",type="contrast");bp1c <- visreg(bp_model,"PPFD_a",type="contrast");
bp1d <- visreg(bp_model,"age_a",type="contrast");bp1e <- visreg(bp_model,"CNrt_a",type="contrast")

w1 <- ggplot(data = bp1a$fit) +geom_line(aes(Tg_a, visregFit),color="black",size=2) + xlab("Tg") + ylab("Measured BP")+
  geom_point(data=bp1a$res,aes(x=Tg_a,y=visregRes),alpha=0.5)+theme_classic()+theme(text = element_text(size=20))+
  geom_ribbon(aes(Tg_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.3)
w2 <- ggplot(data = bp1b$fit) +geom_line(aes(fAPAR_a, visregFit),color="black",size=2) + xlab("fAPAR") + ylab("")+
  geom_point(data=bp1b$res,aes(x=fAPAR_a,y=visregRes),alpha=0.5)+theme_classic()+theme(text = element_text(size=20))+
  geom_ribbon(aes(fAPAR_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.3)
w3 <- ggplot(data = bp1c$fit) +geom_line(aes(PPFD_a, visregFit),color="black",size=2) + xlab("ln PPFD") + ylab("")+
  geom_point(data=bp1c$res,aes(x=PPFD_a,y=visregRes),alpha=0.5)+theme_classic()+theme(text = element_text(size=20))+
  geom_ribbon(aes(PPFD_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.3)
w4 <- ggplot(data = bp1d$fit) +geom_line(aes(age_a, visregFit),color="black",size=2) + xlab("ln age") + ylab("")+
  geom_point(data=bp1d$res,aes(x=age_a,y=visregRes),alpha=0.5)+theme_classic()+theme(text = element_text(size=20))+
  geom_ribbon(aes(age_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.3)
w5 <- ggplot(data = bp1e$fit) +geom_line(aes(CNrt_a, visregFit),color="black",size=2) + xlab("ln soil C/N") + ylab("")+
  geom_point(data=bp1e$res,aes(x=CNrt_a,y=visregRes),alpha=0.5)+theme_classic()+theme(text = element_text(size=20))+
  geom_ribbon(aes(CNrt_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.3)

#trendy model
NPP_statistical <- merge(NPP_forest,sitemean_final,by=c("lon","lat"),all.x=TRUE)
NPP_statistical <- subset(NPP_statistical,TNPP_1>0)

d1 <- na.omit(NPP_statistical[,c("lon","lat","CABLE_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d1<- aggregate(d1,by=list(d1$lon,d1$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t1 <- stepwise_lm(d1,"CABLE_NPP")
t1[[2]]
mod1 <- (lm(CABLE_NPP~fAPAR_a+Tg_a+vpd_a,data=d1))
summary(mod1)

d2 <- na.omit(NPP_statistical[,c("lon","lat","ISAM_npp","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d2<- aggregate(d2,by=list(d2$lon,d2$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t2 <- stepwise_lm(d2,"ISAM_npp")
t2[[1]]
t2[[2]]
mod2 <- (lm(ISAM_npp~fAPAR_a+Tg_a+PPFD_a+vpd_a,data=d2))
summary(mod2)

d3 <- na.omit(NPP_statistical[,c("lon","lat","ISBA_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d3<- aggregate(d3,by=list(d3$lon,d3$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t3 <- stepwise_lm(d3,"ISBA_NPP")
t3[[2]]
mod3 <- (lm(ISBA_NPP~fAPAR_a+Tg_a+vpd_a,data=d3))
summary(mod3)

d4 <- na.omit(NPP_statistical[,c("lon","lat","JULES_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d4<- aggregate(d4,by=list(d4$lon,d4$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t4 <- stepwise_lm(d4,"JULES_NPP")
t4[[2]]
mod4 <- (lm(JULES_NPP~Tg_a+PPFD_a,data=d4))
summary(mod4)

d5 <- na.omit(NPP_statistical[,c("lon","lat","LPJ_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d5<- aggregate(d5,by=list(d5$lon,d5$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t5 <- stepwise_lm(d5,"LPJ_NPP")
t5[[2]]
mod5<- (lm(LPJ_NPP~fAPAR_a+PPFD_a,data=d5))
summary(mod5)

d6 <- na.omit(NPP_statistical[,c("lon","lat","ORCHIDEE_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d6<- aggregate(d6,by=list(d6$lon,d6$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t6 <- stepwise_lm(d6,"ORCHIDEE_NPP")
t6[[2]]
mod6<- (lm(ORCHIDEE_NPP~fAPAR_a+Tg_a+PPFD_a,data=d6))
summary(mod6)

d7 <- na.omit(NPP_statistical[,c("lon","lat","ORCHICNP_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d7<- aggregate(d7,by=list(d7$lon,d7$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t7 <- stepwise_lm(d7,"ORCHICNP_NPP")
t7[[2]]
mod7 <- (lm(ORCHICNP_NPP~fAPAR_a+Tg_a+vpd_a,data=d7))
summary(mod7)

d8 <- na.omit(NPP_statistical[,c("lon","lat","SDGVM_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d8<- aggregate(d8,by=list(d8$lon,d8$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t8 <- stepwise_lm(d8,"SDGVM_NPP")
t8[[2]]
mod8<- (lm(SDGVM_NPP~PPFD_a+Tg_a+vpd_a,data=d8))
summary(mod8)

mod1$coefficients;mod2$coefficients;mod3$coefficients;mod4$coefficients;
mod5$coefficients;mod6$coefficients;mod7$coefficients;mod8$coefficients

t1a <- visreg(mod1,"Tg_a",type="contrast");t2a <-visreg(mod2,"Tg_a",type="contrast");t3a <-visreg(mod3,"Tg_a",type="contrast");t4a <- visreg(mod4,"Tg_a",type="contrast");t6a <-visreg(mod6,"Tg_a",type="contrast");t7a <- visreg(mod7,"Tg_a",type="contrast");t8a <-visreg(mod8,"Tg_a",type="contrast")
f1a <- visreg(mod1,"fAPAR_a",type="contrast");f2a <-visreg(mod2,"fAPAR_a",type="contrast");f3a <-visreg(mod3,"fAPAR_a",type="contrast");f5a <-visreg(mod5,"fAPAR_a",type="contrast");f6a <-visreg(mod6,"fAPAR_a",type="contrast");f7a <- visreg(mod7,"fAPAR_a",type="contrast")
p2a <- visreg(mod2,"PPFD_a",type="contrast");p4a <- visreg(mod4,"PPFD_a",type="contrast");p5a <- visreg(mod5,"PPFD_a",type="contrast");p6a <- visreg(mod6,"PPFD_a",type="contrast");p8a <- visreg(mod8,"PPFD_a",type="contrast")
v1a <- visreg(mod1,"vpd_a",type="contrast");v2a <- visreg(mod2,"vpd_a",type="contrast");v3a <- visreg(mod3,"vpd_a",type="contrast");v7a <- visreg(mod7,"vpd_a",type="contrast");v8a <- visreg(mod8,"vpd_a",type="contrast")

fits_tg <- dplyr::bind_rows(mutate(bp1a$fit, plt = "Measurement"),mutate(t1a$fit, plt = "CABLE"),mutate(t2a$fit, plt = "ISAM"),mutate(t3a$fit, plt = "ISBA"),mutate(t4a$fit, plt = "JULES"),mutate(t6a$fit, plt = "ORCHIDEE"),mutate(t7a$fit, plt = "ORCHICNP"),mutate(t8a$fit, plt = "SDGVM"))

fits_fapar <- dplyr::bind_rows(mutate(bp1b$fit, plt = "Measurement"),mutate(f1a$fit, plt = "CABLE"),mutate(f2a$fit, plt = "ISAM"),mutate(f3a$fit, plt = "ISBA"),mutate(f5a$fit, plt = "LPJ"),mutate(f6a$fit, plt = "ORCHIDEE"),mutate(f7a$fit, plt = "ORCHICNP"))

fits_PPFD <- dplyr::bind_rows(mutate(bp1c$fit, plt = "Measurement"),mutate(p2a$fit, plt = "ISAM"),mutate(p4a$fit, plt = "JULES"),mutate(p5a$fit, plt = "LPJ"),mutate(p6a$fit, plt = "ORCHIDEE"),mutate(p8a$fit, plt = "SDGVM"))

fits_vpd <- dplyr::bind_rows(mutate(v1a$fit, plt = "CABLE"),mutate(v2a$fit, plt = "ISAM"),mutate(v3a$fit, plt = "ISBA"),mutate(v7a$fit, plt = "ORCHICNP"),mutate(v8a$fit, plt = "SDGVM"))

final1 <- ggplot() +geom_line(data = fits_tg, aes(Tg_a, visregFit, group=plt, color=plt),size=2) + xlab("Tg") + ylab("TRENDY BP")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+  geom_ribbon(data = bp1a$fit,aes(Tg_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))
final2 <- ggplot() +geom_line(data = fits_fapar, aes(fAPAR_a, visregFit, group=plt, color=plt),size=2) + xlab("fAPAR") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+geom_ribbon(data = bp1b$fit,aes(fAPAR_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+ scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))
final3 <- ggplot() +geom_line(data = fits_vpd, aes(vpd_a, visregFit, group=plt, color=plt),size=2)+ xlab("ln D") + ylab(" ") +theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))
final4 <- ggplot() +geom_line(data = fits_PPFD, aes(PPFD_a, visregFit, group=plt, color=plt),size=2) + xlab("ln PPFD") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+geom_ribbon(data = bp1c$fit,aes(PPFD_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.5)+  scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))

#show legend
final1_legend <- ggplot() +geom_line(data = fits_PPFD, aes(PPFD_a, visregFit, group=plt, color=plt),size=2) + xlab("ln PPFD") + ylab(" ")+theme_classic()+theme(text = element_text(size=20))+ scale_colour_manual(" ",values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))

legend_info <- as_ggplot(get_legend(final1_legend))

plot_grid(w1,w2,w3,w4,w5,
          final1,final2,final4,final3,legend_info,
          nrow=2,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/newphy_figs3.jpg",sep=""), width = 20, height = 10)

#nuptake figure
#N minerlization
Nmin_statistical <- subset(NPP_forest,Nmin>0)

sitemean2 <- unique(Nmin_statistical[,c("lon","lat")])
sp_sites2 <- SpatialPoints(sitemean2) # only select lon and lat

sitemean2_final <- sitemean2

allmaps2 <- list(ORCHICNP_fNup,ISAM_fNup)
obs_name2 <- c("ORCHICNP_fNup","ISAM_fNup")

for (i in c(1:length(allmaps2))){
  df1 <- allmaps2[[i]]
  names(df1) <- c("lon","lat",obs_name2[i])
  coordinates(df1) <- ~lon+lat 
  gridded(df1) <- TRUE
  df1_global <- raster(df1, obs_name2[i]) 
  
  sp_sites_new <- raster::extract(df1_global, sp_sites2, sp = TRUE) %>% as_tibble() %>% 
    right_join(sitemean2, by = c("lon", "lat"))
  sitemean2_final[,i+2] <- sp_sites_new[,1]
}

sitemean2_final$ORCHICNP_fNup[sitemean2_final$ORCHICNP_fNup==0] <- NA
sitemean2_final$ORCHICNP_fNup[sitemean2_final$ORCHICNP_fNup==0] <- NA

Nmin_statistical <- merge(Nmin_statistical,sitemean2_final,by=c("lon","lat"),all.x=TRUE)
Nmin_statistical$Nmin_a <- Nmin_statistical$Nmin
Nmin_statistical_final <- na.omit(Nmin_statistical[,c("Nmin_a","age_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                                                      "CNrt_a","LMA_a","site_a")])

#Nup model
n1 <- stepwise(Nmin_statistical_final,"Nmin_a")
n1[[1]]
n1[[2]]
summary(lmer(Nmin_a~fAPAR_a+Tg_a+(1|site_a),data=Nmin_statistical_final))
r.squaredGLMM(lmer(Nmin_a~fAPAR_a+Tg_a+(1|site_a),data=Nmin_statistical_final))

#TRENDY1
m1 <- na.omit(Nmin_statistical[,c("lon","lat","ORCHICNP_fNup","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
m1<- aggregate(m1,by=list(m1$lon,m1$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
q1 <- stepwise_lm(m1,"ORCHICNP_fNup")
q1[[2]]
summary(lm(ORCHICNP_fNup~fAPAR_a,data=m1)) #remove vpd as it only consist of Tg and fapar

m2 <- na.omit(Nmin_statistical[,c("lon","lat","ISAM_fNup","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
m2<- aggregate(m2,by=list(m2$lon,m2$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
q2 <- stepwise_lm(m2,"ISAM_fNup")
q2[[2]]
summary(lm(ISAM_fNup~fAPAR_a+vpd_a,data=m2))

#final figure for Nup
mod_n1 <- lmer(Nmin_a~fAPAR_a+Tg_a+(1|site_a),data=Nmin_statistical_final)
r.squaredGLMM(mod_n1)
summary(mod_n1)
nn1a <- visreg(mod_n1,"fAPAR_a",type="contrast")
nn1b <- visreg(mod_n1,"Tg_a",type="contrast")

ww1 <- ggplot(data = nn1a$fit) +geom_line(aes(fAPAR_a, visregFit),color="black",size=2) + xlab("fAPAR") + ylab("Measured Nmin")+
  geom_point(data=nn1a$res,aes(x=fAPAR_a,y=visregRes),alpha=0.5)+theme_classic()+theme(text = element_text(size=20))+
  geom_ribbon(aes(fAPAR_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.3)

ww2 <- ggplot(data = nn1b$fit) +geom_line(aes(Tg_a, visregFit),color="black",size=2) + xlab("Tg") + ylab("Measured Nmin")+
  geom_point(data=nn1b$res,aes(x=Tg_a,y=visregRes),alpha=0.5)+theme_classic()+theme(text = element_text(size=20))+
  geom_ribbon(aes(Tg_a, ymin=visregLwr, ymax=visregUpr),fill="gray",alpha=0.3)

mod_n2 <- lm(ORCHICNP_fNup~fAPAR_a,data=m1)
nn1c <- visreg(mod_n2,"fAPAR_a",type="contrast")

mod_n3 <- lm(ISAM_fNup~fAPAR_a,data=m2)
nn1d <- visreg(mod_n3,"fAPAR_a",type="contrast")

fits_fapar <- dplyr::bind_rows(mutate(nn1a$fit, plt = "Measurement"),mutate(nn1c$fit, plt = "ORCHICNP"),mutate(nn1d$fit, plt = "ISAM_fNup"))
#fits_tg <- dplyr::bind_rows(mutate(nn1b$fit, plt = "Measurement"),mutate(nn1d$fit, plt = "ISAM"))

final1a <- ggplot() +geom_line(data = fits_fapar, aes(fAPAR_a, visregFit, group=plt, color=plt),size=2) + xlab("fAPAR") + ylab("TRENDY N uptake")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(Measurement="black",ISAM="#663399",ORCHICNP="cyan"))
#final1b <- ggplot() +geom_line(data = fits_tg, aes(Tg_a, visregFit, group=plt, color=plt),size=2) + xlab("Tg") + ylab("TRENDY N uptake")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(Measurement="black",ISAM="#663399",ORCHICNP="cyan"))

final2_legend <- ggplot() +geom_line(data = fits_fapar, aes(fAPAR_a, visregFit, group=plt, color=plt),size=2) + xlab("fAPAR") + ylab("TRENDY N uptake")+theme_classic()+theme(text = element_text(size=20))+ scale_colour_manual(" ",values=c(Measurement="black",ISAM="#663399",ORCHICNP="cyan"))

legend_info2 <- as_ggplot(get_legend(final2_legend))

plot_grid(ww1,ww2,white,
          final1a,legend_info2,
          nrow=2,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/newphy_figs4.jpg",sep=""), width = 15, height = 10)

#nue model for TRENDY
#nuptake figure
#N minerlization
NUE_test <- subset(NPP_all,file!="Tiandi Grassland")

sitemean3 <- unique(NUE_test[,c("lon","lat")])
sp_sites3 <- SpatialPoints(sitemean3) # only select lon and lat

sitemean3_final <- sitemean3

allmaps3 <- list(ORCHICNP_fNup,ISAM_fNup,ORCHICNP_NPP,ISAM_npp)
obs_name3 <- c("ORCHICNP_fNup","ISAM_fNup","ORCHICNP_NPP","ISAM_npp")

for (i in c(1:length(allmaps3))){
  df1 <- allmaps3[[i]]
  names(df1) <- c("lon","lat",obs_name3[i])
  coordinates(df1) <- ~lon+lat 
  gridded(df1) <- TRUE
  df1_global <- raster(df1, obs_name3[i]) 
  
  sp_sites_new <- raster::extract(df1_global, sp_sites3, sp = TRUE) %>% as_tibble() %>% 
    right_join(sitemean3, by = c("lon", "lat"))
  sitemean3_final[,i+2] <- sp_sites_new[,1]
}

sitemean3_final$ORCHICNP_fNup[sitemean3_final$ORCHICNP_fNup==0] <- NA
sitemean3_final$ORCHICNP_fNup[sitemean3_final$ORCHICNP_fNup==0] <- NA
sitemean3_final$ISAM_npp[sitemean3_final$ISAM_npp==0] <- NA

NUE_test <- merge(NUE_test,sitemean3_final,by=c("lon","lat"),all.x=TRUE)
NUE_test$NUE_ORCHICNP <- NUE_test$ORCHICNP_NPP/NUE_test$ORCHICNP_fNup
NUE_test$NUE_ISAM <- NUE_test$ISAM_npp/NUE_test$ISAM_fNup

NUE_test_sitemean <- aggregate(NUE_test,by=list(NUE_test$lon,NUE_test$lat), FUN=mean, na.rm=TRUE)

NUE_test_final <- na.omit(NUE_test_sitemean[,c("NUE_ORCHICNP","age_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                                                      "CNrt_a","LMA_a","vcmax25_a")])

NUE_test_final2 <- na.omit(NUE_test_sitemean[,c("NUE_ISAM",
                                      "age_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                                      "CNrt_a","LMA_a","vcmax25_a")])


#NUE_ORCHICNP
nue1 <- stepwise_lm(NUE_test_final,"NUE_ORCHICNP")
nue1[[1]]
nue1[[2]]
summary(lm(NUE_ORCHICNP~LMA_a+fAPAR_a+age_a+vpd_a+CNrt_a+Tg_a+vcmax25_a,data=NUE_test_final))

#NUE_ISAM
nue2 <- stepwise_lm(NUE_test_final2,"NUE_ISAM")
nue2[[1]]
nue2[[2]]
summary(lm(NUE_ISAM~CNrt_a+Tg_a+age_a+vpd_a+PPFD_a,data=NUE_test_final2))


#validation - BP
NPP_statistical$Measured_BP <- NPP_statistical$TNPP_1

pp2 <- analyse_modobs2(NPP_statistical,"pred_npp","Measured_BP", type = "points")
pp4 <- analyse_modobs2(NPP_statistical,"CABLE_NPP","Measured_BP", type = "points")
pp5 <- analyse_modobs2(NPP_statistical,"ISAM_npp","Measured_BP", type = "points")
pp6 <- analyse_modobs2(NPP_statistical,"ISBA_NPP","Measured_BP", type = "points")
pp7 <- analyse_modobs2(NPP_statistical,"JULES_NPP","Measured_BP", type = "points")
pp8 <- analyse_modobs2(NPP_statistical,"LPJ_NPP","Measured_BP", type = "points")
pp9 <- analyse_modobs2(NPP_statistical,"ORCHIDEE_NPP","Measured_BP", type = "points")
pp10 <- analyse_modobs2(NPP_statistical,"ORCHICNP_NPP","Measured_BP", type = "points")
pp11 <- analyse_modobs2(NPP_statistical,"SDGVM_NPP","Measured_BP", type = "points")

plot_grid(pp2$gg,pp4$gg,pp5$gg,pp6$gg,pp7$gg,pp8$gg,pp9$gg,pp10$gg,pp11$gg,
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'),
          nrow=3,label_x = 0.9,label_y=0.92)+white

ggsave(paste("~/data/output/newphy_figs5.jpg",sep=""), width = 15, height = 16)

#validation - N uptake
Nmin_statistical$pred_nuptake

Nmin_statistical$Nuptake_from_Nup_model <- summary(mod_n1)$coef[1,1]+ summary(mod_n1)$coef[2,1]*Nmin_statistical$fAPAR_a+
  summary(mod_n1)$coef[3,1]*Nmin_statistical$Tg_a

ppp1 <- analyse_modobs2(Nmin_statistical,"pred_nuptake","Nmin", type = "points")
ppp2 <- analyse_modobs2(Nmin_statistical,"Nuptake_from_Nup_model","Nmin", type = "points")
ppp3 <- analyse_modobs2(Nmin_statistical,"ORCHICNP_fNup","Nmin", type = "points")
ppp4 <- analyse_modobs2(Nmin_statistical,"ISAM_fNup","Nmin", type = "points")

plot_grid(ppp1$gg,ppp2$gg,ppp3$gg,ppp4$gg,
          labels = c('(a)','(b)','(c)','(d)'),
          nrow=2,label_x = 0.9,label_y=0.92)+white
ggsave(paste("~/data/output/newphy_figs6.jpg",sep=""), width = 10, height = 10)

#visreg

a1 <- ~{
  p1a <- visreg(bp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a2 <- ~{
  p1a <- visreg(bp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a3 <- ~{
  p1a <- visreg(bp_model,"CNrt_a",type="contrast")
  plot(p1a,ylab="Forest BP",xlab="ln soil C/N",
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
  plot(p1a,ylab="logit ANPP/BP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a7 <- ~{
  p1a <- visreg(anpp_tnpp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a8 <- ~{
  p1a <- visreg(anpp_tnpp_model,"CNrt_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln soil C/N",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a9 <- ~{
  p1a <- visreg(anpp_tnpp_model,"age_a",type="contrast")
  plot(p1a,ylab="logit ANPP/BP",xlab="ln age",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a10 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"Tg_a",type="contrast")
  plot(p1a,ylab="logit leaf-NPP/ANPP",xlab="Tg",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a11 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"PPFD_a",type="contrast")
  plot(p1a,ylab="logit leaf-NPP/ANPP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}

a12 <- ~{
  p1a <- visreg(anpp_leafnpp_model,"vpd_a",type="contrast")
  plot(p1a,ylab="logit leaf-NPP/ANPP",xlab="ln D",
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
  plot(p1a,ylab="Grassland BP",xlab="ln PPFD",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)}


plot_grid(a1,a2,a3,a4,a5,
          a6,a7,a8,a9,white,
          a10,a11,a12,white,white,
          a13,a14,white,white,white,
          a15,a16,white,white,white,
          nrow=5)+white

ggsave(paste("~/data/output/newphy_figs1.jpg",sep=""), width = 20, height = 20)

#table s1
names(all_maps)
for (i in 4:31){
  varname <- names(all_maps)[i]
  total_value <- round(sum(all_maps[,i]*conversion,na.rm=TRUE),2)
  print(varname)
  print(total_value)
}
