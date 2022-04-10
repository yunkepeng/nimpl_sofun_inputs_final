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
rm(list=ls())
#stepwise function, where site_a is a random factor
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
#1. trendy model output


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
plot_map3(na.omit(JULES_GPP[,c("lon","lat","myvar")]),varnam = "myvar",latmin = -65, latmax = 85)

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
sitemean <- unique(read.csv("~/data/NPP_final/NPP_validation.csv")[,c("lon","lat")])
sp_sites <- SpatialPoints(sitemean) # only select lon and lat

sitemean_final <- unique(read.csv("~/data/NPP_final/NPP_validation.csv")[,c("lon","lat")])

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

#
summary(sitemean_final) #convert 0 to NA
sitemean_final$CABLE_NPP[sitemean_final$CABLE_NPP==0] <- NA
sitemean_final$ISAM_fNup[sitemean_final$ISAM_fNup==0] <- NA
sitemean_final$ISAM_gpp[sitemean_final$ISAM_gpp==0] <- NA
sitemean_final$ISAM_npp[sitemean_final$ISAM_npp==0] <- NA
sitemean_final$JULES_GPP[sitemean_final$JULES_GPP==0] <- NA
sitemean_final$JULES_NPP[sitemean_final$JULES_NPP==0] <- NA
summary(sitemean_final) 

validation <- read.csv("~/data/NPP_final/NPP_validation.csv")
NPP_statistical <- merge(validation,sitemean_final,by=c("lon","lat"),all.x=TRUE)

NPP_statistical$obs_age[NPP_statistical$obs_age==999] <- NA

NPP_statistical$tnpp_a <- NPP_statistical$TNPP_1

NPP_statistical$soilCN_a <- log(NPP_statistical$soilCN)
NPP_statistical$observedfAPAR_a <- NPP_statistical$observedfAPAR
NPP_statistical$obs_age_a <- log(NPP_statistical$obs_age)

NPP_statistical$age_a <- log(NPP_statistical$age)
NPP_statistical$alpha_a <- (NPP_statistical$alpha)
NPP_statistical$Tg_a <- NPP_statistical$Tg
NPP_statistical$PPFD_a <- log(NPP_statistical$PPFD)
NPP_statistical$vpd_a <- log(NPP_statistical$vpd)
NPP_statistical$fAPAR_a <- NPP_statistical$fAPAR
NPP_statistical$CNrt_a <- log(NPP_statistical$CNrt)
NPP_statistical$LMA_a <- log(NPP_statistical$LMA)
NPP_statistical$vcmax25_a <- log(NPP_statistical$max_vcmax25_c3)
NPP_statistical$site_a <- NPP_statistical$site

#examine alternative model using mapped predictor - too bad!
a1 <- stepwise(na.omit(NPP_statistical[,c("tnpp_gpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")]),"tnpp_gpp_a")
a1[[1]]

a1_old <- stepwise(na.omit(NPP_statistical[,c("tnpp_gpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")]),"tnpp_gpp_a")
a1_old[[1]]
summary(lmer(tnpp_gpp_a~obs_age_a+observedfAPAR_a+soilCN_a+(1|site_a),data=NPP_statistical))

a2 <- stepwise(na.omit(NPP_statistical[,c("anpp_tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")]),"anpp_tnpp_a")
a2[[1]]

r.squaredGLMM(lmer(tnpp_gpp_a~CNrt_a+(1|site_a),data=NPP_statistical))
r.squaredGLMM(lmer(anpp_tnpp_a~CNrt_a+(1|site_a),data=NPP_statistical))

#two dataset
#1. large - all mapping data
NPP_statistical_large <- na.omit(NPP_statistical[,c("tnpp_a","age_a","Tg_a","PPFD_a","vpd_a","fAPAR_a","CNrt_a","site_a")])

#2. small - observed data
NPP_statistical_small <- na.omit(NPP_statistical[,c("tnpp_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])

# r2 =0.44 at large dataset 
a1<- stepwise(NPP_statistical_large,"tnpp_a")
a1[[1]]
model_1 <-  (lmer(tnpp_a~Tg_a+fAPAR_a+age_a+CNrt_a+PPFD_a+(1|site_a),data=NPP_statistical_large))
r.squaredGLMM(model_1)

# r2 =0.68 at small 
a3<- stepwise(NPP_statistical_small,"tnpp_a")
a3[[3]]
model_2 <- (lmer(tnpp_a~Tg_a+observedfAPAR_a+obs_age_a+PPFD_a+(1|site_a),data=NPP_statistical_small))
r.squaredGLMM(model_2)

#1st additional figure 
white <- theme(plot.background=element_rect(fill="white", color="white"))

bp1a <- visreg(model_1,"Tg_a",type="contrast");bp1b <- visreg(model_1,"fAPAR_a",type="contrast");bp1c <- visreg(model_1,"PPFD_a",type="contrast");
bp1d <- visreg(model_1,"age_a",type="contrast");bp1e <- visreg(model_1,"CNrt_a",type="contrast")

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


####1. Large dataset
#trendy's larger dataset - using lm! Because response BNPP are at site-level. It also has higher R2
#Large dataset age_a,alpha_a,Tg_a,PPFD_a,vpd_a,fAPAR_a,CNrt_a

#only available measured npp points - shall we?
NPP_statistical <- subset(NPP_statistical,TNPP_1>0)

d1 <- na.omit(NPP_statistical[,c("lon","lat","CABLE_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d1<- aggregate(d1,by=list(d1$lon,d1$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t1 <- stepwise_lm(d1,"CABLE_NPP")
t1[[2]]
mod1 <- (lm(CABLE_NPP~fAPAR_a+Tg_a+vpd_a+PPFD_a,data=d1))

d2 <- na.omit(NPP_statistical[,c("lon","lat","ISAM_npp","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d2<- aggregate(d2,by=list(d2$lon,d2$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t2 <- stepwise_lm(d2,"ISAM_npp")
t2[[2]]
mod2 <- (lm(ISAM_npp~fAPAR_a+Tg_a+PPFD_a,data=d2))

d3 <- na.omit(NPP_statistical[,c("lon","lat","ISBA_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d3<- aggregate(d3,by=list(d3$lon,d3$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t3 <- stepwise_lm(d3,"ISBA_NPP")
t3[[2]]
mod3 <- (lm(ISBA_NPP~fAPAR_a+Tg_a+vpd_a,data=d3))

d4 <- na.omit(NPP_statistical[,c("lon","lat","JULES_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d4<- aggregate(d4,by=list(d4$lon,d4$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t4 <- stepwise_lm(d4,"JULES_NPP")
t4[[2]]
mod4 <- (lm(JULES_NPP~Tg_a+PPFD_a,data=d4))

d5 <- na.omit(NPP_statistical[,c("lon","lat","LPJ_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d5<- aggregate(d5,by=list(d5$lon,d5$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t5 <- stepwise_lm(d5,"LPJ_NPP")
t5[[2]]
mod5<- (lm(LPJ_NPP~fAPAR_a+PPFD_a,data=d5))

d6 <- na.omit(NPP_statistical[,c("lon","lat","ORCHIDEE_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d6<- aggregate(d6,by=list(d6$lon,d6$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t6 <- stepwise_lm(d6,"ORCHIDEE_NPP")
t6[[2]]
mod6<- (lm(ORCHIDEE_NPP~fAPAR_a+Tg_a+PPFD_a,data=d6))

d7 <- na.omit(NPP_statistical[,c("lon","lat","ORCHICNP_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d7<- aggregate(d7,by=list(d7$lon,d7$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t7 <- stepwise_lm(d7,"ORCHICNP_NPP")
t7[[2]]
mod7 <- (lm(ORCHICNP_NPP~fAPAR_a+Tg_a+vpd_a,data=d7))

d8 <- na.omit(NPP_statistical[,c("lon","lat","SDGVM_NPP","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
d8<- aggregate(d8,by=list(d8$lon,d8$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t8 <- stepwise_lm(d8,"SDGVM_NPP")
t8[[2]]
mod8<- (lm(SDGVM_NPP~PPFD_a+Tg_a+vpd_a,data=d8))

####2. Small dataset
d1a <- na.omit(NPP_statistical[,c("lon","lat","CABLE_NPP","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d1a<- aggregate(d1a,by=list(d1a$lon,d1a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t1a <- stepwise_lm(d1a,"CABLE_NPP")
t1a[[2]]
summary(lm(CABLE_NPP~soilCN_a+Tg_a+PPFD_a,data=d1a))

d2a <- na.omit(NPP_statistical[,c("lon","lat","ISAM_npp","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d2a<- aggregate(d2a,by=list(d2a$lon,d2a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t2a <- stepwise_lm(d2a,"ISAM_npp")
t2a[[2]]
summary(lm(ISAM_npp~Tg_a+PPFD_a+observedfAPAR_a,data=d2a))

d3a <- na.omit(NPP_statistical[,c("lon","lat","ISBA_NPP","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d3a<- aggregate(d3a,by=list(d3a$lon,d3a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t3a <- stepwise_lm(d3a,"ISBA_NPP")
t3a[[2]]
summary(lm(ISBA_NPP~Tg_a+alpha_a+observedfAPAR_a,data=d3a))

d4a <- na.omit(NPP_statistical[,c("lon","lat","JULES_NPP","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d4a <- aggregate(d4a,by=list(d4a$lon,d4a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t4a <- stepwise_lm(d4a,"JULES_NPP")
t4a[[2]]
summary(lm(JULES_NPP~Tg_a+observedfAPAR_a+soilCN_a,data=d4a))

d5a <- na.omit(NPP_statistical[,c("lon","lat","LPJ_NPP","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d5a <- aggregate(d5a,by=list(d5a$lon,d5a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t5a <- stepwise_lm(d5a,"LPJ_NPP")
t5a[[2]]
summary(lm(LPJ_NPP~Tg_a+vpd_a+soilCN_a,data=d5a))

d6a <- na.omit(NPP_statistical[,c("lon","lat","ORCHIDEE_NPP","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d6a<- aggregate(d6a,by=list(d6a$lon,d6a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t6a <- stepwise_lm(d6a,"ORCHIDEE_NPP")
t6a[[2]]
summary(lm(ORCHIDEE_NPP~Tg_a+PPFD_a+observedfAPAR_a+vpd_a,data=d6a))

d7a <- na.omit(NPP_statistical[,c("lon","lat","ORCHICNP_NPP","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d7a <- aggregate(d7a,by=list(d7a$lon,d7a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t7a <- stepwise_lm(d7a,"ORCHICNP_NPP")
t7a[[2]]
summary(lm(ORCHICNP_NPP~PPFD_a+observedfAPAR_a+alpha_a+Tg_a+soilCN_a,data=d7a))

d8a <- na.omit(NPP_statistical[,c("lon","lat","SDGVM_NPP","alpha_a","Tg_a","PPFD_a","vpd_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")])
d8a <- aggregate(d8a,by=list(d8a$lon,d8a$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
t8a <- stepwise_lm(d8a,"SDGVM_NPP")
t8a[[2]]
summary(lm(SDGVM_NPP~Tg_a+PPFD_a+vpd_a,data=d8a))



#3rd additional figure - fig
#model  trendy BP fitted by mapped predictors
mod1$coefficients;mod2$coefficients;mod3$coefficients;mod4$coefficients;
mod5$coefficients;mod6$coefficients;mod7$coefficients;mod8$coefficients

t1a <- visreg(mod1,"Tg_a",type="contrast");t2a <-visreg(mod2,"Tg_a",type="contrast");t3a <-visreg(mod3,"Tg_a",type="contrast");t4a <- visreg(mod4,"Tg_a",type="contrast");t6a <-visreg(mod6,"Tg_a",type="contrast");t7a <- visreg(mod7,"Tg_a",type="contrast");t8a <-visreg(mod8,"Tg_a",type="contrast")
f1a <- visreg(mod1,"fAPAR_a",type="contrast");f2a <-visreg(mod2,"fAPAR_a",type="contrast");f3a <-visreg(mod3,"fAPAR_a",type="contrast");f5a <-visreg(mod5,"fAPAR_a",type="contrast");f6a <-visreg(mod6,"fAPAR_a",type="contrast");f7a <- visreg(mod7,"fAPAR_a",type="contrast")
p1a <- visreg(mod1,"PPFD_a",type="contrast");p2a <- visreg(mod2,"PPFD_a",type="contrast");p4a <- visreg(mod4,"PPFD_a",type="contrast");p5a <- visreg(mod5,"PPFD_a",type="contrast");p6a <- visreg(mod6,"PPFD_a",type="contrast");p8a <- visreg(mod8,"PPFD_a",type="contrast")
v1a <- visreg(mod1,"vpd_a",type="contrast");v3a <- visreg(mod3,"vpd_a",type="contrast");v7a <- visreg(mod7,"vpd_a",type="contrast");v8a <- visreg(mod8,"vpd_a",type="contrast")

fits_tg <- dplyr::bind_rows(mutate(bp1a$fit, plt = "Measurement"),mutate(t1a$fit, plt = "CABLE"),mutate(t2a$fit, plt = "ISAM"),mutate(t3a$fit, plt = "ISBA"),mutate(t4a$fit, plt = "JULES"),mutate(t6a$fit, plt = "ORCHIDEE"),mutate(t7a$fit, plt = "ORCHICNP"),mutate(t8a$fit, plt = "SDGVM"))

fits_fapar <- dplyr::bind_rows(mutate(bp1b$fit, plt = "Measurement"),mutate(f1a$fit, plt = "CABLE"),mutate(f2a$fit, plt = "ISAM"),mutate(f3a$fit, plt = "ISBA"),mutate(f5a$fit, plt = "LPJ"),mutate(f6a$fit, plt = "ORCHIDEE"),mutate(f7a$fit, plt = "ORCHICNP"))

fits_PPFD <- dplyr::bind_rows(mutate(bp1c$fit, plt = "Measurement"),mutate(p1a$fit, plt = "CABLE"),mutate(p2a$fit, plt = "ISAM"),mutate(p4a$fit, plt = "JULES"),mutate(p5a$fit, plt = "LPJ"),mutate(p8a$fit, plt = "SDGVM"))

fits_vpd <- dplyr::bind_rows(mutate(v1a$fit, plt = "CABLE"),mutate(v3a$fit, plt = "ISBA"),mutate(v7a$fit, plt = "ORCHICNP"),mutate(v8a$fit, plt = "SDGVM"))

final1 <- ggplot() +geom_line(data = fits_tg, aes(Tg_a, visregFit, group=plt, color=plt),size=2) + xlab("Tg") + ylab("TRENDY BP")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))
final2 <- ggplot() +geom_line(data = fits_fapar, aes(fAPAR_a, visregFit, group=plt, color=plt),size=2) + xlab("fAPAR") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))
final3 <- ggplot() +geom_line(data = fits_vpd, aes(vpd_a, visregFit, group=plt, color=plt),size=2)+ xlab("ln D") + ylab(" ") +theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))
final4 <- ggplot() +geom_line(data = fits_PPFD, aes(PPFD_a, visregFit, group=plt, color=plt),size=2) + xlab("ln PPFD") + ylab(" ")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))

#show legend
final1_legend <- ggplot() +geom_line(data = fits_PPFD, aes(PPFD_a, visregFit, group=plt, color=plt),size=2) + xlab("ln PPFD") + ylab(" ")+theme_classic()+theme(text = element_text(size=20))+ scale_colour_manual(" ",values=c(Measurement="black",CABLE="pink",ISAM="#663399",ISBA="#339999",JULES="#CC0033",LPJ="#FF6600",ORCHIDEE="#FF9933",ORCHICNP="cyan",SDGVM="yellow"))

legend_info <- as_ggplot(get_legend(final1_legend))

plot_grid(w1,w2,w3,w4,w5,
          final1,final2,final4,final3,legend_info,
          nrow=2,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/figs_new1.jpg",sep=""), width = 20, height = 10)


#N minerlization
Nmin_statistical <- read.csv("/Users/yunpeng/data/NPP_final/Nmin_validation.csv")
Nmin_statistical <- subset(Nmin_statistical,Nmin>0)

sitemean2 <- unique(read.csv("/Users/yunpeng/data/NPP_final/Nmin_validation.csv")[,c("lon","lat")])
sp_sites2 <- SpatialPoints(sitemean2) # only select lon and lat

sitemean2_final <- unique(read.csv("/Users/yunpeng/data/NPP_final/Nmin_validation.csv")[,c("lon","lat")])

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

Nmin_statistical$age_a <- log(Nmin_statistical$age)
Nmin_statistical$alpha_a <- (Nmin_statistical$alpha)
Nmin_statistical$Tg_a <- Nmin_statistical$Tg
Nmin_statistical$PPFD_a <- log(Nmin_statistical$PPFD)
Nmin_statistical$vpd_a <- log(Nmin_statistical$vpd)
Nmin_statistical$fAPAR_a <- Nmin_statistical$fAPAR
Nmin_statistical$CNrt_a <- log(Nmin_statistical$CNrt)
Nmin_statistical$LMA_a <- log(Nmin_statistical$LMA)
Nmin_statistical$vcmax25_a <- log(Nmin_statistical$max_vcmax25_c3)
Nmin_statistical$site_a <- Nmin_statistical$sitename

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
summary(lm(ORCHICNP_fNup~fAPAR_a,data=m1))

m2 <- na.omit(Nmin_statistical[,c("lon","lat","ISAM_fNup","Tg_a","PPFD_a","vpd_a","fAPAR_a","site_a")])
m2<- aggregate(m2,by=list(m2$lon,m2$lat), FUN=mean, na.rm=TRUE)%>% dplyr::select(-c(Group.1,Group.2,lon,lat,site_a))
q2 <- stepwise_lm(m2,"ISAM_fNup")
q2[[2]]
summary(lm(ISAM_fNup~Tg_a,data=m2))

#final figure for Nup
mod_n1 <- lmer(Nmin_a~fAPAR_a+Tg_a+(1|site_a),data=Nmin_statistical_final)

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

mod_n3 <- lm(ISAM_fNup~Tg_a,data=m2)
nn1d <- visreg(mod_n3,"Tg_a",type="contrast")

fits_fapar <- dplyr::bind_rows(mutate(nn1a$fit, plt = "Measurement"),mutate(nn1c$fit, plt = "ORCHICNP"))
fits_tg <- dplyr::bind_rows(mutate(nn1b$fit, plt = "Measurement"),mutate(nn1d$fit, plt = "ISAM"))

final1a <- ggplot() +geom_line(data = fits_fapar, aes(fAPAR_a, visregFit, group=plt, color=plt),size=2) + xlab("fAPAR") + ylab("TRENDY N uptake")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(Measurement="black",ISAM="#663399",ORCHICNP="cyan"))
final1b <- ggplot() +geom_line(data = fits_tg, aes(Tg_a, visregFit, group=plt, color=plt),size=2) + xlab("Tg") + ylab("TRENDY N uptake")+theme_classic()+theme(text = element_text(size=20),legend.position="none")+ scale_colour_manual(values=c(Measurement="black",ISAM="#663399",ORCHICNP="cyan"))

final2_legend <- ggplot() +geom_line(data = fits_tg, aes(Tg_a, visregFit, group=plt, color=plt),size=2) + xlab("Tg") + ylab("TRENDY N uptake")+theme_classic()+theme(text = element_text(size=20))+ scale_colour_manual(" ",values=c(Measurement="black",ISAM="#663399",ORCHICNP="cyan"))

legend_info2 <- as_ggplot(get_legend(final2_legend))

plot_grid(ww1,ww2,white,
          final1a,final1b,legend_info2,
          nrow=2,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/figs_new2.jpg",sep=""), width = 15, height = 10)

#validation - BP
NPP_statistical$pred_npp_model1 <- 
  summary(model_1)$coef[1,1]+summary(model_1)$coef[2,1]*NPP_statistical$Tg_a+summary(model_1)$coef[3,1]*NPP_statistical$fAPAR_a+
  summary(model_1)$coef[4,1]*NPP_statistical$age_a+summary(model_1)$coef[5,1]*NPP_statistical$CNrt_a+summary(model_1)$coef[6,1]*NPP_statistical$PPFD_a

NPP_statistical$BP_from_BPE_model <- NPP_statistical$pred_npp
NPP_statistical$BP_from_BP_model <- NPP_statistical$pred_npp_model1

NPP_statistical$Measured_BP <- NPP_statistical$TNPP_1

p1 <- analyse_modobs2(NPP_statistical,"BP_from_BPE_model","Measured_BP", type = "points")
p2 <- analyse_modobs2(NPP_statistical,"BP_from_BP_model","Measured_BP", type = "points")
p4 <- analyse_modobs2(NPP_statistical,"CABLE_NPP","Measured_BP", type = "points")
p5 <- analyse_modobs2(NPP_statistical,"ISAM_npp","Measured_BP", type = "points")
p6 <- analyse_modobs2(NPP_statistical,"ISBA_NPP","Measured_BP", type = "points")
p7 <- analyse_modobs2(NPP_statistical,"JULES_NPP","Measured_BP", type = "points")
p8 <- analyse_modobs2(NPP_statistical,"LPJ_NPP","Measured_BP", type = "points")
p9 <- analyse_modobs2(NPP_statistical,"ORCHIDEE_NPP","Measured_BP", type = "points")
p10 <- analyse_modobs2(NPP_statistical,"ORCHICNP_NPP","Measured_BP", type = "points")
p11 <- analyse_modobs2(NPP_statistical,"SDGVM_NPP","Measured_BP", type = "points")

plot_grid(p1$gg,p2$gg,p4$gg,p5$gg,p6$gg,p7$gg,p8$gg,p9$gg,p10$gg,p11$gg,
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)'),
          nrow=4,label_x = 0.9,label_y=0.92)+white

ggsave(paste("~/data/output/figs_new3.jpg",sep=""), width = 15, height = 16)

#bpe
NPP_statistical$BPE_from_BPE_model <- NPP_statistical$pred_npp/NPP_statistical$pred_gpp_c3
NPP_statistical$BPE_from_BP_model <- NPP_statistical$pred_npp_model1/NPP_statistical$pred_gpp_c3

NPP_statistical$CABLE_BPE <- NPP_statistical$CABLE_NPP/NPP_statistical$CABLE_GPP
NPP_statistical$ISAM_BPE <- NPP_statistical$ISAM_npp/NPP_statistical$ISAM_gpp
NPP_statistical$ISBA_BPE <- NPP_statistical$ISBA_NPP/NPP_statistical$ISBA_GPP
NPP_statistical$JULES_BPE <- NPP_statistical$JULES_NPP/NPP_statistical$JULES_GPP
NPP_statistical$LPJ_BPE <- NPP_statistical$LPJ_NPP/NPP_statistical$LPJ_GPP
NPP_statistical$ORCHIDEE_BPE <- NPP_statistical$ORCHIDEE_NPP/NPP_statistical$ORCHIDEE_GPP
NPP_statistical$ORCHICNP_BPE <- NPP_statistical$ORCHICNP_NPP/NPP_statistical$ORCHICNP_GPP
NPP_statistical$SDGVM_BPE <- NPP_statistical$SDGVM_NPP/NPP_statistical$SDGVM_GPP
NPP_statistical$Measured_BPE <- NPP_statistical$TNPP_1/NPP_statistical$GPP

p1a <- analyse_modobs2(NPP_statistical,"BPE_from_BPE_model","Measured_BPE", type = "points")
#p2a <- analyse_modobs2(NPP_statistical,"BPE_from_BP_model","Measured_BPE", type = "points")
p4a <- analyse_modobs2(NPP_statistical,"CABLE_BPE","Measured_BPE", type = "points")
p5a <- analyse_modobs2(NPP_statistical,"ISAM_BPE","Measured_BPE", type = "points")
p6a <- analyse_modobs2(NPP_statistical,"ISBA_BPE","Measured_BPE", type = "points")
p7a <- analyse_modobs2(NPP_statistical,"JULES_BPE","Measured_BPE", type = "points")
p8a <- analyse_modobs2(NPP_statistical,"LPJ_BPE","Measured_BPE", type = "points")
p9a <- analyse_modobs2(NPP_statistical,"ORCHIDEE_BPE","Measured_BPE", type = "points")
p10a <- analyse_modobs2(NPP_statistical,"ORCHICNP_BPE","Measured_BPE", type = "points")
p11a <- analyse_modobs2(NPP_statistical,"SDGVM_BPE","Measured_BPE", type = "points")

plot_grid(p1a$gg,p4a$gg,p5a$gg,p6a$gg,p7a$gg,p8a$gg,p9a$gg,p10a$gg,p11a$gg,
          labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'),
          nrow=3,label_x = 0.9,label_y=0.92)+white

ggsave(paste("~/data/output/figs_new4.jpg",sep=""), width = 15, height = 16)

#validation - N uptake
Nmin_statistical$Nuptake_from_MS <-Nmin_statistical$pred_nuptake
summary(mod_n1)$coef

Nmin_statistical$Nuptake_from_Nup_model <- summary(mod_n1)$coef[1,1]+ summary(mod_n1)$coef[2,1]*Nmin_statistical$fAPAR_a+
  summary(mod_n1)$coef[3,1]*Nmin_statistical$Tg_a

pp1 <- analyse_modobs2(Nmin_statistical,"Nuptake_from_MS","Nmin", type = "points")
pp2 <- analyse_modobs2(Nmin_statistical,"Nuptake_from_Nup_model","Nmin", type = "points")
pp3 <- analyse_modobs2(Nmin_statistical,"ORCHICNP_fNup","Nmin", type = "points")
pp4 <- analyse_modobs2(Nmin_statistical,"ISAM_fNup","Nmin", type = "points")

plot_grid(pp1$gg,pp2$gg,pp3$gg,pp4$gg,
          labels = c('(a)','(b)','(c)','(d)'),
          nrow=2,label_x = 0.9,label_y=0.92)+white

ggsave(paste("~/data/output/figs_new5.jpg",sep=""), width = 10, height = 10)


#global!
####copied from global_figs.R
#add existing files
###1. load all prediction fields map (details in ~/yunkepeng/nimpl_sofun_inputs_final/Prediction_field), gpp and vcmax25 (derived from SOFUN/yunkebranch).
firstyr_data <- 1982 # In data file, which is the first year
endyr_data <- 2011 # In data file, which is the last year
location <- "~/data/output/latest_forest/"
alloutput_list <- list.files(location,full.names = T)

#input elevation nc file, which will be cbind with global df directly
elev_nc <- read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc")
#elev_nc <- read_nc_onefile("D:/PhD/nimpl_sofun_inputs/Data/Elevation/WFDEI-elevation.nc")
elev <- as.data.frame(nc_to_df(elev_nc, varnam = "elevation"))

#2. Create a function to specify path, loop many years nc file and output a dataframe (lon, lat, var).
inputnc <- function(name,start_year,end_year){
  #-----------------------------------------------------------------------
  # Input: 
  # name: gpp, npp, anpp, vcmax25, leafcn, nuptake...
  # start_year: e.g. 1981
  # end_year: e.g. 2016
  # location: e.g "D:/PhD/nimpl_sofun_inputs/Data/output/" or in Euler: "~/yunkebranch_units/outputnc/"
  #-----------------------------------------------------------------------
  output_allyears <- data.frame(matrix(NA))
  # first, include all years annual data into a daframe
  for (i in firstyr_data:endyr_data){
    if (name == "npp"){
      nc <- read_nc_onefile(alloutput_list[grepl("a.npp.nc", list.files(location,full.names = T))][i-firstyr_data+1]) #we only rely this to filter npp.nc file...
    } else {
      nc <- read_nc_onefile(alloutput_list[grepl(name, list.files(location,full.names = T))][i-firstyr_data+1]) #Input nc
    }
    output_year <- nc_to_df(nc, varnam = name)[,3] #Yearly output
    output_allyears[1:259200,i-firstyr_data+1] <- output_year #here first column represents first year of data file 's output
  }
  names(output_allyears) <- paste(name,firstyr_data:endyr_data,sep="")
  #this variable above (output_allyears), could be end of the function, which is variable at multiple years. But for our purporses, we need mean of select years
  #then, only calculate means of selected years
  output_selected_yrs <- rowMeans(output_allyears[,(start_year-firstyr_data+1):(end_year-firstyr_data+1)],na.rm = TRUE) # only calculated means based on selected start and end year (see function)
  coord <- nc_to_df(nc, varnam = name)[,1:2] # obtain lon and lat
  final_output <- cbind(coord,elev[,3],output_selected_yrs) # combine lon, lat,z with rowmeans variable
  names(final_output) <- c("lon","lat","z",name)
  return(final_output)
  #-----------------------------------------------------------------------
  # Output: output_final: the output data (259200 * 3) including lon, lat and value
  #-----------------------------------------------------------------------
}

#3. select data over 30 years, each df includes lon, lat, z, var
vcmax25_df <- inputnc("vcmax25",1982,2011)

gpp_df <- inputnc("gpp",1982,2011)

#now, inputting all predictors
Tg <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/Tg.nc"),
  varnam = "Tg"))

PPFD <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/PPFD.nc"),
  varnam = "PPFD"))

vpd <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/vpd.nc"),
  varnam = "vpd"))
#alpha not used...
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

nre_constant_grass <- 0.69

#input all regressions 

#firstly, load all forest models
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

#now, do the same for grassland
load(file = "~/data/NPP_Grassland_final/statistical_model/tnpp_grass.RData")
mod_tnpp_grass<- tnpp_grass
summary(mod_tnpp_grass)
load(file = "~/data/NPP_Grassland_final/statistical_model/anpp_grass.RData")
mod_anpp_grass <- anpp_grass
summary(mod_anpp_grass)

###2. run global simulations
npp_f <- gpp_df$gpp * (1/(1 + exp(-(summary(mod_tnpp)$coef[1,1]+
                                      summary(mod_tnpp)$coef[2,1]* log(CNrt$myvar)+
                                      summary(mod_tnpp)$coef[3,1] * log(age$myvar) + 
                                      summary(mod_tnpp)$coef[4,1]* fAPAR$myvar))))

anpp_f <- gpp_df$gpp * (1/(1 + exp(-(summary(mod_anpp)$coef[1,1]+
                                       summary(mod_anpp)$coef[2,1] * log(CNrt$myvar)+ 
                                       summary(mod_anpp)$coef[3,1] * log(age$myvar) + 
                                       summary(mod_anpp)$coef[4,1] * fAPAR$myvar))))

bnpp_f <- npp_f-anpp_f

lnpp_f <- anpp_f * (1/(1 + exp(-(summary(mod_lnpp)$coef[1,1]+
                                   summary(mod_lnpp)$coef[2,1]* log(PPFD$myvar) +
                                   summary(mod_lnpp)$coef[3,1] * (Tg$myvar) +
                                   summary(mod_lnpp)$coef[4,1] * log(vpd$myvar)))))

wnpp_f <- anpp_f - lnpp_f

leafnc_f <- (summary(n1)$coef[1,1]/0.46) + 
  (summary(n1)$coef[2,1]/0.46) *vcmax25_df$vcmax25/LMA$myvar
#0.46 is constant Cmass

nre_f <- (1/(1+exp(-(summary(nre_model)$coef[1,1]+
                       summary(nre_model)$coef[2,1] *Tg$myvar + 
                       summary(nre_model)$coef[3,1] * log(vpd$myvar)))))

lnf_f <- (1-nre_f)* leafnc_f * lnpp_f

wnf_f <- wnpp_f/100
#100 is constant wood c/n

bnf_f <- bnpp_f/94
#94 is constant root c/n

nuptake_f <- lnf_f + wnf_f + bnf_f

#grass
npp_g <- gpp_df$gpp * summary(mod_tnpp_grass)$coef[1,1]
anpp_g <- gpp_df$gpp * summary(mod_anpp_grass)$coef[1,1]
bnpp_g <- npp_g-anpp_g

leafnc_g <- 1/18

nre_g <- 0.69

lnf_g <- anpp_g*leafnc_g*(1-nre_g)

bnf_g <- bnpp_g *(1/41)
#41 is constant root c/n
nuptake_g <- lnf_g + bnf_g

###2. input land cover
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
                                      fAPAR$myvar,age$myvar,
                                      CNrt$myvar,LMA$myvar,vcmax25_df$vcmax25))
all_predictors$available_grid = rowMeans(all_predictors)
#just to find all na columns
all_predictors$available_grid[is.na(all_predictors$available_grid)==FALSE] <- 1
summary(all_predictors$available_grid)
available_grid2 <- all_predictors$available_grid

#represent grids when stand-age is especially in NA, but others are fine
names(all_predictors) <- c("Tg","PPFD","vpd","fAPAR","age","CNrt","LMA","vcmax25","available_grid")
all_predictors$lon <- gpp_df$lon
all_predictors$lat <- gpp_df$lat
summary(all_predictors)


#final calculation - now divide into forest, grassland and pft
#available_grid2 here was used as a list of data to identify if a grid is available (=1) or any prediction fields shown as NA 
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


all_maps <- as.data.frame(cbind(gpp_df,npp_pft,npp_forest,npp_grass,
                                anpp_pft,anpp_forest,anpp_grass,
                                bnpp_pft,bnpp_forest,bnpp_grass,
                                lnpp_forest,wnpp_forest,wnf_forest,
                                leafcn_pft,leafcn_forest,leafcn_grassland,
                                nre_pft,nre_forest,nre_grassland,
                                lnf_pft,lnf_forest,lnf_grass,
                                bnf_pft,bnf_forest,bnf_grass,
                                nuptake_pft,nuptake_forest,nuptake_grass))

summary(all_maps)

all_maps$CUE <- all_maps$npp_pft/all_maps$gpp
all_maps$NUE <- all_maps$npp_pft/all_maps$nuptake_pft

b5 <- ggplot(data=all_maps, aes(x=CUE, y=NUE)) +
  geom_point(aes(color=nuptake_pft),alpha=0.3,size=0.3)+geom_smooth(method = "lm", se = TRUE)+ 
  scale_color_viridis(discrete=FALSE,direction= -1)+theme_classic()+theme(axis.title = element_text(size = 20),
                                                                          axis.text = element_text(size = 15),
                                                                          legend.title = element_text(size = 14))+
  xlab("BPE")+ylab("NUE (gC/gN)")+labs(color= ~paste("N uptake", " (gN m"^-2,"yr"^-1,")"))

summary(lm(all_maps$NUE~all_maps$CUE))

#now, trendy
ORCHICNP_final <- as.data.frame(cbind(ORCHICNP_fNup[,3],ORCHICNP_GPP[,3],ORCHICNP_NPP[,3]))
names(ORCHICNP_final) <- c("Nuptake","gpp","npp")
ISAM_final <-as.data.frame(cbind(ISAM_fNup[,3],ISAM_gpp[,3],ISAM_npp[,3]))
names(ISAM_final) <- c("Nuptake","gpp","npp")

ORCHICNP_final$gpp[ORCHICNP_final$gpp==0] <-NA
ORCHICNP_final$npp[ORCHICNP_final$npp==0] <-NA
ORCHICNP_final$Nuptake[ORCHICNP_final$Nuptake==0] <-NA
ORCHICNP_final$CUE <- ORCHICNP_final$npp/ORCHICNP_final$gpp
ORCHICNP_final$NUE <- ORCHICNP_final$npp/ORCHICNP_final$Nuptake

ISAM_final$gpp[ISAM_final$gpp==0] <-NA
ISAM_final$npp[ISAM_final$npp==0] <-NA
ISAM_final$Nuptake[ISAM_final$Nuptake==0] <-NA
ISAM_final$CUE <- ISAM_final$npp/ISAM_final$gpp
ISAM_final$NUE <- ISAM_final$npp/ISAM_final$Nuptake

b6 <- ggplot(data=ORCHICNP_final, aes(x=CUE, y=NUE)) +xlim(c(0,1))+ylim(c(0,200))+
  geom_point(aes(color=Nuptake),alpha=0.3,size=0.3)+geom_smooth(method = "lm", se = TRUE)+ 
  scale_color_viridis(discrete=FALSE,direction= -1)+theme_classic()+theme(axis.title = element_text(size = 20),
                                                                          axis.text = element_text(size = 15),
                                                                          legend.title = element_text(size = 14))+
  xlab("ORCHICNP BPE")+ylab("ORCHICNP NUE (gC/gN)")+labs(color= ~paste("N uptake", " (gN m"^-2,"yr"^-1,")"))

summary(lm(ISAM_final$NUE~ISAM_final$CUE))

b7 <- ggplot(data=ISAM_final, aes(x=CUE, y=NUE)) +xlim(c(0,1))+ylim(c(0,200))+
  geom_point(aes(color=Nuptake),alpha=0.3,size=0.3)+geom_smooth(method = "lm", se = TRUE)+ 
  scale_color_viridis(discrete=FALSE,direction= -1)+theme_classic()+theme(axis.title = element_text(size = 20),
                                                                          axis.text = element_text(size = 15),
                                                                          legend.title = element_text(size = 14))+
  xlab("ISAM BPE")+ylab("ISAM NUE (gC/gN)")+labs(color= ~paste("N uptake", " (gN m"^-2,"yr"^-1,")"))

summary(lm(ISAM_final$NUE~ISAM_final$CUE))

plot_grid(b5,b6,b7, nrow=1,label_x = 0.8, label_y = 0.8)+white

ggsave(paste("~/data/output/figs_new6.jpg",sep=""), width = 20, height = 5)
