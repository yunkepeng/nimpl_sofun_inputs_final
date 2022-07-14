#check siteinfo

#check overall
NPP_all <- read.csv("~/data/NPP_Yunke/NPP_Nmin_dataset_with_predictors.csv")
dim(unique(NPP_all[,c("lon","lat","z")]))
length(unique(NPP_all[,c("site")]))
check <- (unique(NPP_all[,c("lon","lat","z","site")]))

(unique(NPP_all[,c("file")]))

subset_data <- subset(NPP_all,file=="MCampioli")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))
#all unique - sitename number is equal to total numbers
#e.g. US-bol-D01, US-bol-D02 is two different sites


subset_data <- subset(NPP_all,file=="ForC")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))

subset_data <- subset(NPP_all,file=="Sara Vicca")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))

subset_data <- subset(NPP_all,file=="Malhi 2011")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))

#finzi ignored

# Vicca_validation_file - no changed - as it is created by myself
subset_data <- subset(NPP_all,file=="Vicca_validation_file")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))

# Keith  - naerly no change
subset_data <- subset(NPP_all,file=="Keith")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))

# Malhi - created by myself
subset_data <- subset(NPP_all,file=="Malhi 2017")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))

# tiandi - created by myself
subset_data <- subset(NPP_all,file=="~/data/npp_stoichiometry_forests_tiandi/")
dim(subset_data)
nrow(unique(subset_data[,c("lon","lat","z")]))
length(unique(subset_data[,c("site")]))
aa <- (unique(subset_data[,c("lon","lat","z","site")]))

#summary
#try if recreating original name
#BP model didn't change
#anpp/bp model changed to, age effect was removed now
#leaf-npp/anpp changed to C/N (-) + age(+)

library(MLmetrics)
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
library(car)
library("ggplotify")
library(remotes)
library(tune)
library(relaimpo)

#reset validation metrics info

analyse_modobs2 <- function(
    df,
    mod,
    obs,
    type       = "points",
    filnam     = NA,
    relative   = FALSE,
    xlim       = NULL,
    ylim       = NULL,
    use_factor = NULL,
    shortsubtitle = FALSE,
    plot_subtitle = TRUE,
    plot_linmod = TRUE,
    ...
){
  
  require(ggplot2)
  require(dplyr)
  require(LSD)
  require(ggthemes)
  require(RColorBrewer)
  
  #if (identical(filnam, NA)) filnam <- "analyse_modobs.pdf"
  
  ## rename to 'mod' and 'obs' and remove rows with NA in mod or obs
  df <- df %>%
    as_tibble() %>%
    ungroup() %>%
    dplyr::select(mod=mod, obs=obs) %>%
    tidyr::drop_na(mod, obs)
  
  ## get linear regression (coefficients)
  linmod <- lm( obs ~ mod, data=df )
  
  ## construct metrics table using the 'yardstick' library
  df_metrics <- df %>%
    yardstick::metrics(obs, mod) %>%
    dplyr::bind_rows( tibble( .metric = "n",        .estimator = "standard", .estimate = summarise(df, numb=n()) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "slope",    .estimator = "standard", .estimate = coef(linmod)[2]) ) %>%
    # dplyr::bind_rows( tibble( .metric = "nse",      .estimator = "standard", .estimate = hydroGOF::NSE( obs, mod, na.rm=TRUE ) ) ) %>%
    dplyr::bind_rows( tibble( .metric = "mean_obs", .estimator = "standard", .estimate = summarise(df, mean=mean(obs, na.rm=TRUE)) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "prmse",    .estimator = "standard",
                              .estimate = dplyr::filter(., .metric=="rmse") %>% dplyr::select(.estimate) %>% unlist() /
                                dplyr::filter(., .metric=="mean_obs") %>% dplyr::select(.estimate) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "pmae",    .estimator = "standard",
                              .estimate = dplyr::filter(., .metric=="mae") %>% dplyr::select(.estimate) %>% unlist() /
                                dplyr::filter(., .metric=="mean_obs") %>% dplyr::select(.estimate) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "bias",        .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod-obs), na.rm=TRUE    )) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "pbias",       .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod-obs)/obs, na.rm=TRUE)) %>% unlist() ) )
  
  rsq_val <- df_metrics %>% dplyr::filter(.metric=="rsq") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  rmse_val <- df_metrics %>% dplyr::filter(.metric=="rmse") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  mae_val <- df_metrics %>% dplyr::filter(.metric=="mae") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  bias_val <- df_metrics %>% dplyr::filter(.metric=="bias") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  slope_val <- df_metrics %>% dplyr::filter(.metric=="slope") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  n_val <- df_metrics %>% dplyr::filter(.metric=="n") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  
  if (relative){
    rmse_val <- rmse_val / mean(df$obs, na.rm = TRUE)
    bias_val <- bias_val / mean(df$obs, na.rm = TRUE)
  }
  
  rsq_lab <- format( rsq_val, digits = 2 )
  rmse_lab <- format( rmse_val, digits = 3 )
  mae_lab <- format( mae_val, digits = 3 )
  bias_lab <- format( bias_val, digits = 3 )
  slope_lab <- format( slope_val, digits = 3 )
  n_lab <- format( n_val, digits = 3 )
  
  results <- tibble( rsq = rsq_val, rmse = rmse_val, mae = mae_val, bias = bias_val, slope = slope_val, n = n_val )
  
  if (shortsubtitle){
    subtitle <- bquote( italic(R)^2 == .(rsq_lab) ~~
                          RRMSE == .(rmse_lab) )
  } else {
    subtitle <- bquote( italic(R)^2 == .(rsq_lab) ~~
                          RRMSE == .(rmse_lab) ~~
                          bias == .(bias_lab) ~~
                          slope == .(slope_lab) ~~
                          italic(N) == .(n_lab) )
  }
  
  if (type=="heat"){
    
    # if (!identical(filnam, NA)) dev.off()
    # source("~/LSD/R/LSD.heatscatter.R")
    
    gg <- heatscatter(
      df$mod,
      df$obs,
      xlim=xlim,
      ylim=ylim,
      main="",
      ggplot=TRUE )
    
    gg <- gg +
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="hex"){
    
    ## ggplot hexbin
    gg <- df %>%
      ggplot2::ggplot(aes(x=mod, y=obs)) +
      geom_hex() +
      scale_fill_gradientn(
        colours = colorRampPalette( c("gray65", "navy", "red", "yellow"))(5)) +
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="points") {
    
    ## points
    gg <- df %>%
      ggplot(aes(x=mod, y=obs)) +
      geom_point() +
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="density") {
    
    ## points
    gg <- df %>%
      ggplot(aes(x=mod, y=obs)) +
      
      stat_density_2d(aes(fill = after_stat(nlevel)), geom = "polygon") +
      scale_fill_gradientn(colours = colorRampPalette( c("gray65", "navy", "red", "yellow"))(5),
                           guide = "legend") +
      
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  }
  
  return(list(df_metrics=df_metrics, gg=gg, linmod=linmod, results = results))
}

white <- theme(plot.background=element_rect(fill="white", color="white"))

larger_size <- theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                     plot.subtitle=element_text(size=15))

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

#re-create site-name
aa <- unique(NPP_all[,c("lon","lat","z")])
aa$sitename <- paste("a",1:nrow(aa),sep="")

#merge into site
NPP_all <- merge(NPP_all,aa,by=c("lon","lat","z"),all.x=TRUE)

NPP_all$site <- NPP_all$sitename

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

#check why some grassland BP and ANPP is so high
outliers <- subset(NPP_all,pft=="Grassland" & ANPP_2>500)
newmap <- getMap(resolution = "low")
sp::plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(outliers$lon,outliers$lat, col="green", pch=16,cex=2)

BP_dataset <- na.omit(NPP_forest[,c("tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")])
#model1 <- stepwise(BP_dataset,"tnpp_a")
#model1[[1]]
#model1[[2]]
#bp_model <- (lmer(tnpp_a~Tg_a+observedfAPAR_a+obs_age_a+PPFD_a+alpha_a+(1|site_a),data=BP_dataset))
#summary(bp_model)

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

#anpp_tnpp_dataset <- na.omit(NPP_forest[,c("anpp_tnpp_a","obs_age_a","observedfAPAR_a","soilCN_a","Tg_a","PPFD_a","vpd_a","site_a")])
#dim(anpp_tnpp_dataset)
#model2 <- stepwise(anpp_tnpp_dataset,"anpp_tnpp_a")
#model2[[1]]
#model2[[2]]
#anpp_tnpp_model <- (lmer(anpp_tnpp_a~soilCN_a+obs_age_a+observedfAPAR_a+(1|site_a),data=anpp_tnpp_dataset))
#summary(anpp_tnpp_model)

#mapped
anpp_tnpp_dataset2 <- na.omit(NPP_forest[,c("anpp_tnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
dim(subset(NPP_forest,ANPP_2>0))
model2a <- stepwise(anpp_tnpp_dataset2,"anpp_tnpp_a")
model2a[[1]]
model2a[[2]]
model2a[[3]]
anpp_tnpp_model <- (lmer(anpp_tnpp_a~CNrt_a+PPFD_a+Tg_a+(1|site_a),data=anpp_tnpp_dataset2))
summary(anpp_tnpp_model)
r.squaredGLMM(anpp_tnpp_model)

vif_anpp_tnpp <- vif((lmer(anpp_tnpp_a~Tg_a+fAPAR_a+PPFD_a+CNrt_a+age_a+vpd_a+(1|site_a),data=anpp_tnpp_dataset2)))

#with age, but age should be removed since it shows higher AIC and lower R2
anpp_leafnpp_dataset_age <- na.omit(NPP_forest[,c("anpp_leafnpp_a","age_a","fAPAR_a","CNrt_a","Tg_a","PPFD_a","vpd_a","site_a")])
model3a <- stepwise(anpp_leafnpp_dataset_age,"anpp_leafnpp_a")
model3a[[1]]
model3a[[2]]
model3a[[3]]
test <- (lmer(anpp_leafnpp_a~Tg_a+(1|site_a),data=anpp_leafnpp_dataset_age)) 
summary(test)

#without Tg - re-selection - this is best
anpp_leafnpp_dataset <- na.omit(NPP_forest[,c("anpp_leafnpp_a","fAPAR_a","CNrt_a","age_a","PPFD_a","vpd_a","site_a")])
model3 <- stepwise(anpp_leafnpp_dataset,"anpp_leafnpp_a")
model3[[1]]
model3[[3]]
anpp_leafnpp_model <- (lmer(anpp_leafnpp_a~CNrt_a+age_a+(1|site_a),data=anpp_leafnpp_dataset)) 
summary(anpp_leafnpp_model)
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
model_g1[[3]]

bp_grass_model <- (lm(tnpp_a~PPFD_a+Tg_a,data=BP_dataset_grass))
summary(bp_grass_model)
r.squaredGLMM(bp_grass_model)

vif_bp_grass <- vif((lm(tnpp_a~Tg_a+PPFD_a+vpd_a+CNrt_a+fAPAR_a,data=BP_dataset_grass)))

#anpp/tnpp
dim(subset(grassland_sitemean,TNPP_1>0))
dim(subset(grassland_sitemean,TNPP_1>0 & ANPP_2>0))

anpp_tnpp_dataset_grass <- na.omit(grassland_sitemean[,c("anpp_tnpp_a","Tg_a","PPFD_a","vpd_a","CNrt_a","fAPAR_a")])
model_g2 <- stepwise_lm(anpp_tnpp_dataset_grass,"anpp_tnpp_a")
model_g2[[1]]
summary(lm(anpp_tnpp_a~Tg_a,data=anpp_tnpp_dataset_grass))
#non-significant! so alternatively using constant ratio
summary((lm(anpp_a~-1+tnpp_a,data=grassland_sitemean))) # 0.49 for anpp, so 0.51 for bnpp
