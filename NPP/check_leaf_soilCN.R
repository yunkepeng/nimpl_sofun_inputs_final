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

#modify function
#' Analyse modelled values versus observed data.
#'
#' Calculates a set of performance statistics and optionally creates plots of modelled
#' versus observed values.
#'
#' @param df A data frame containing columns with names corresponding to arguments
#' \code{mod} and \code{obs}
#' @param mod A character string specifying the variable name (column) of the
#' modelled (simulated) values in data frame \code{df}.
#' @param obs A character string specifying the variable name (column) of the
#' observed values in data frame \code{df}.
#' @param type If \code{"points"}, uses \code{geom_points()}, if \code{"hex"}
#' uses \code{ggplot2::geom_hex()}, if \code{"heat"} uses adjusted
#' \code{geom_points()} with color indicating density, if \code{"density"} uses
#' \code{stat_density_2d()} to draws polygos of equal density.
#' @param filnam A character string specifying the name of the file containing
#' the plot. Defaults to \code{NA} (no file is created).
#' @param relative A logical specifying whether the relative RMSE and bias (after
#' division by the mean) is to be showed in the subtitle labels.
#' @param shortsubtitle A boolean specifying whether to display a reduced set of metrics
#' in the subtitle.
#' @param rsquared A boolean specifying whether to display R-squared and the RMSE
#' (if \code{TRUE}) or the r (Pearson's correlation coefficient) and the p (p-value of
#' test of significance of correlation, if \code{TRUE}). Defaluts to \code{TRUE}.
#' @param plot_subtitle A boolean specifying whether to display any metrics. Defaults
#' to \code{TRUE}.
#' @param plot_linmod A boolean specifying whether to display the fitted linear
#' regression as a red line. Defaults to \code{TRUE}.
#' @param plot_legend A boolean specifying whether to display a legend for the colors.
#' Defaults to \code{TRUE} if \code{type} is one of  \code{"heat"},  \code{"hex"}, or
#' \code{"density"}.
#' @param label A boolean specifying whether points should be labelled using ggrepel.
#' Defaults to \code{FALSE}. Only available for \code{type == "points"}. Use argument
#' \code{nlabels} to specify how many points should be labelled, starting with points
#' that have the largest residuals from the linear regression fit.
#' @param id A character string specifying the column name that identifies the points.
#' The column's values must be of type integer and is used to label points in case of
#' \code{label = TRUE}.
#' @param nlabels An integer specifying how many points to be labelled, starting with points
#' that have the largest residuals from the linear regression fit. Only available
#' for \code{type == "points"}. Defaults to one.
#'
#' @export
#'
#' @examples
#'
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
    rsquared    = TRUE,
    plot_subtitle = TRUE,
    plot_linmod = TRUE,
    plot_legend = TRUE,
    label       = FALSE,
    id          = NULL,
    nlabels     = 1,
    ...
){
  
  require(ggplot2)
  require(dplyr)
  require(LSD)
  require(ggthemes)
  require(RColorBrewer)
  
  #if (identical(filnam, NA)) filnam <- "analyse_modobs.pdf"
  
  ## rename to 'mod' and 'obs' and remove rows with NA in mod or obs
  if (label){
    df <- df %>%
      as_tibble() %>%
      ungroup() %>%
      dplyr::select(mod=mod, obs=obs, id=!!id) %>%
      tidyr::drop_na(mod, obs)
    
  } else {
    df <- df %>%
      as_tibble() %>%
      ungroup() %>%
      dplyr::select(mod=mod, obs=obs) %>%
      tidyr::drop_na(mod, obs)
    
  }
  
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
    dplyr::bind_rows( tibble( .metric = "bias",  .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod-obs), na.rm=TRUE    )) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "pbias", .estimator = "standard", .estimate = dplyr::summarise(df, mean((mod-obs)/obs, na.rm=TRUE)) %>% unlist() ) ) %>%
    dplyr::bind_rows( tibble( .metric = "cor",   .estimator = "standard", .estimate = cor(df$mod, df$obs, method = "pearson") ) ) %>%
    dplyr::bind_rows( tibble( .metric = "cor_p", .estimator = "standard", .estimate = cor.test(df$mod, df$obs, method = "pearson")$p.value ) )
  
  rsq_val <- df_metrics %>% dplyr::filter(.metric=="rsq") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  rmse_val <- df_metrics %>% dplyr::filter(.metric=="rmse") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  mae_val <- df_metrics %>% dplyr::filter(.metric=="mae") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  bias_val <- df_metrics %>% dplyr::filter(.metric=="bias") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  slope_val <- df_metrics %>% dplyr::filter(.metric=="slope") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  n_val <- df_metrics %>% dplyr::filter(.metric=="n") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  cor_val <- df_metrics %>% dplyr::filter(.metric=="cor") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  cor_p_val <- df_metrics %>% dplyr::filter(.metric=="cor_p") %>% dplyr::select(.estimate) %>% unlist() %>% unname()
  
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
  cor_lab <- format( cor_val, digits = 3 )
  cor_p_lab <- format( cor_p_val, digits = 3 )
  
  results <- tibble( rsq = rsq_val, rmse = rmse_val, mae = mae_val, bias = bias_val, slope = slope_val, n = n_val )
  
  if (shortsubtitle){
    if (rsquared){
      subtitle <- bquote(
        italic(R)^2 == .(rsq_lab) ~~
          RMSE == .(rmse_lab)
      )
    } else {
      subtitle <- bquote(
        italic(r) == .(cor_lab) ~~
          italic(p) == .(cor_p_lab)
      )
    }
  } else {
    subtitle <- bquote( italic(R)^2 == .(rsq_lab) ~~
                          RMSE == .(rmse_lab) ~~
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
    
    if (plot_linmod) gg <- gg #+ geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (!plot_legend) gg <- gg + theme(legend.position = "none")
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="hex"){
    
    ## ggplot hexbin
    gg <- df %>%
      ggplot2::ggplot(aes(x=mod, y=obs)) +
      geom_hex(bins = 100) +
      scale_fill_gradientn(
        colours = colorRampPalette( c("gray65", "navy", "red", "yellow"))(5),
        trans = "log") +
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg #+ geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    if (!plot_legend) gg <- gg + theme(legend.position = "none")
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="points"){
    
    if (label){
      df <- df %>%
        dplyr::mutate(.res = mod - obs) %>%
        dplyr::mutate(.absres = abs(.res)) %>%
        dplyr::arrange(desc(.absres)) %>%
        dplyr::mutate(.rankres = 1:n()) %>%
        dplyr::mutate(.dolab = ifelse(.rankres <= nlabels, TRUE, FALSE))
      
      ## points with labels
      library(ggrepel)
      gg <- df %>%
        ggplot(aes(x=mod, y=obs, label = ifelse(.dolab, id, ""))) +
        geom_point() +
        geom_label_repel(min.segment.length = 0, seed = 42, box.padding = 0.5) +
        geom_point(data = dplyr::filter(df, .dolab), color = "red") + #geom_abline(intercept=0, slope=1, linetype="dotted") +
        theme_classic() +
        labs(x = mod, y = obs)
      
    } else {
      ## points
      gg <- df %>%
        ggplot(aes(x=mod, y=obs)) +
        geom_point() +#geom_abline(intercept=0, slope=1, linetype="dotted") +
        theme_classic() +
        labs(x = mod, y = obs)
      
    }
    
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg #+ geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  } else if (type=="density"){
    
    ## density as raster
    gg <- df %>%
      ggplot(aes(x=mod, y=obs)) +
      
      stat_density_2d(
        geom = "raster", #the geometric object to display the data (in this case: rectangles)
        aes(fill = after_stat(density)), #using `density`, a variable calculated by the stat
        contour = FALSE
      ) +
      
      scale_fill_gradientn(colours = colorRampPalette( c("white", "gray65", "navy", "red", "yellow"))(6),
                           guide = FALSE) +
      
      geom_abline(intercept=0, slope=1, linetype="dotted") +
      # coord_fixed() +
      # xlim(0,NA) +
      # ylim(0,NA) +
      theme_classic() +
      labs(x = mod, y = obs)
    
    if (plot_subtitle) gg <- gg + labs(subtitle = subtitle)
    if (plot_linmod) gg <- gg + geom_smooth(method='lm', color="red", size=0.5, se=FALSE)
    if (!plot_legend) gg <- gg + theme(legend.position = "none")
    
    if (!identical(filnam, NA)) {
      ggsave(filnam, width=5, height=5)
    }
    
  }
  
  return(list(df_metrics=df_metrics, gg=gg, linmod=linmod, results = results))
}


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

#output to data frame
csvfile <- paste("~/data/leaf_traits/leaf_soil.csv")
write.csv(SP_input4, csvfile, row.names = TRUE)

SP_input4_sitemean <- aggregate(SP_input4,by=list(SP_input4$lon_2,SP_input4$lat_2), FUN=mean, na.rm=TRUE) #site-mean

#log-transformmed
SP_input4_sitemean$log_narea <- log(SP_input4_sitemean$narea)
SP_input4_sitemean$log_nmass <- log(SP_input4_sitemean$nmass)
SP_input4_sitemean$log_leafCN <- log(SP_input4_sitemean$leafCN)
SP_input4_sitemean$log_LMA <- log(SP_input4_sitemean$lma)
SP_input4_sitemean$log_Vcmax25 <- log(SP_input4_sitemean$Vcmax25)
SP_input4_sitemean$log_cmass <- log(SP_input4_sitemean$cmass)
SP_input4_sitemean$log_carea <- log(SP_input4_sitemean$carea)

larger_size <- theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
                     plot.subtitle=element_text(size=15))

a1 <- analyse_modobs2(SP_input4_sitemean,"CN","log_narea", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil C/N")+larger_size

a2 <- analyse_modobs2(SP_input4_sitemean,"CN","log_nmass", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)

a3 <- analyse_modobs2(SP_input4_sitemean,"CN","log_leafCN", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)

a4 <- analyse_modobs2(SP_input4_sitemean,"CN","log_LMA", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)

a5 <- analyse_modobs2(SP_input4_sitemean,"CN","log_Vcmax25", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil C/N")+larger_size

a6 <- analyse_modobs2(SP_input4_sitemean,"CN","log_cmass", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)

a7 <- analyse_modobs2(SP_input4_sitemean,"CN","log_carea", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil C/N")+larger_size+geom_smooth(method="lm",color="red",size=2)

b1 <- analyse_modobs2(SP_input4_sitemean,"P_total","log_narea", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil P")+larger_size

b2 <- analyse_modobs2(SP_input4_sitemean,"P_total","log_nmass", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil P")+larger_size

b3 <- analyse_modobs2(SP_input4_sitemean,"P_total","log_leafCN", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil P")+larger_size

b4 <- analyse_modobs2(SP_input4_sitemean,"P_total","log_LMA", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil P")+larger_size

b5 <- analyse_modobs2(SP_input4_sitemean,"P_total","log_Vcmax25", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil P")+larger_size

b6 <- analyse_modobs2(SP_input4_sitemean,"P_total","log_cmass", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil P")+larger_size

b7 <- analyse_modobs2(SP_input4_sitemean,"P_total","log_carea", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="Soil P")+larger_size

c1 <- analyse_modobs2(SP_input4_sitemean,"pH","log_narea", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)

c2 <- analyse_modobs2(SP_input4_sitemean,"pH","log_nmass", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2,se=F,linetype = "dashed")

c3 <- analyse_modobs2(SP_input4_sitemean,"pH","log_leafCN", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2,se=F,linetype = "dashed")

c4 <- analyse_modobs2(SP_input4_sitemean,"pH","log_LMA", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)

c5 <- analyse_modobs2(SP_input4_sitemean,"pH","log_Vcmax25", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)

c6 <- analyse_modobs2(SP_input4_sitemean,"pH","log_cmass", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)

c7 <- analyse_modobs2(SP_input4_sitemean,"pH","log_carea", type = "points",shortsubtitle=TRUE,rsquared=FALSE)$gg+labs(x ="pH")+larger_size+geom_smooth(method="lm",color="red",size=2)

plot_grid(a1,a2,a3,a4,a5,
          c1,c2,c3,c4,c5,
          nrow=2)
ggsave(paste("~/data/soilN_P.jpg",sep=""),width = 25,height = 10)

v1 <- (lm(log_cmass~pH+CN,data=SP_input4_sitemean))
summary(v1)

v1a <- visreg(v1,"pH",type="contrast")
v1b <- visreg(v1,"CN",type="contrast")

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
