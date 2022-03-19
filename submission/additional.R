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
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")
#library(rbeni)
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)
NPP_statistical <- read.csv("~/data/NPP_final/NPP_validation.csv")
NPP_statistical$obs_age[NPP_statistical$obs_age==999] <- NA
NPP_statistical <- subset(NPP_statistical,TNPP_1>0)

NPP_statistical$tnpp_a <- NPP_statistical$TNPP_1

NPP_statistical$soilCN_a <- log(NPP_statistical$soilCN)
NPP_statistical$observedfAPAR_a <- NPP_statistical$observedfAPAR
NPP_statistical$obs_age_a <- NPP_statistical$obs_age

NPP_statistical$age_a <- log(NPP_statistical$age)
NPP_statistical$alpha_a <- log(NPP_statistical$alpha)
NPP_statistical$Tg_a <- NPP_statistical$Tg
NPP_statistical$PPFD_a <- log(NPP_statistical$PPFD)
NPP_statistical$vpd_a <- log(NPP_statistical$vpd)
NPP_statistical$fAPAR_a <- NPP_statistical$fAPAR
NPP_statistical$CNrt_a <- log(NPP_statistical$CNrt)
NPP_statistical$LMA_a <- log(NPP_statistical$LMA)
NPP_statistical$vcmax25_a <- log(NPP_statistical$max_vcmax25_c3)
NPP_statistical$site_a <- NPP_statistical$site

#if all -->suggesting npp ~ vcmax25 only is best
NPP_statistical_large <- NPP_statistical[,c("tnpp_a","age_a","alpha_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                                            "CNrt_a","LMA_a","site_a")]

NPP_statistical_middle <- NPP_statistical[,c("tnpp_a","alpha_a","Tg_a","PPFD_a","vpd_a",
                                            "LMA_a","site_a")]

NPP_statistical_small <- NPP_statistical[,c("tnpp_a","alpha_a","Tg_a","PPFD_a","vpd_a",
                                             "LMA_a","site_a","soilCN_a","observedfAPAR_a","obs_age_a")]

NPP_statistical_final <- na.omit(NPP_statistical_small)

# r2 =0.44 at large dataset 
summary(lmer(tnpp_a~Tg_a+fAPAR_a+age_a+CNrt_a+PPFD_a+(1|site_a),data=NPP_statistical_large))

# r2 =0.34 at middle 
summary(lmer(tnpp_a~PPFD_a+Tg_a+(1|site_a),data=NPP_statistical_middle))

# r2 =0.69 at small 
summary(lmer(tnpp_a~Tg_a+observedfAPAR_a+LMA_a+alpha_a+PPFD_a+obs_age_a+soilCN_a+(1|site_a),data=NPP_statistical_small))



#################################
#1. Now, prepare Stepwise regression for NPP/GPP
#################################
library(tidyverse)
library(ggplot2)

#1. select the most important predictor 
summary(NPP_statistical_final)
#determine targets and preds.
target <- 'tnpp_a'

preds <- NPP_statistical_final %>% select(-c(tnpp_a,site_a)) %>% 
  names()

r_list <- c()
library(lme4)
#For loop functions, include all predictor's r2 at the end
for (var in preds){
  forml <- paste( 'lmer(', target, '~', var, '+(1|site_a), data = NPP_statistical_final)')
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

ggplot(All_rsquare, aes(x = reorder(preds, -rsq), y = rsq)) + geom_point() + theme(axis.text.x = element_text(angle = 60, hjust = 1))

#2. stepwise regression selection
library(caret)
library(recipes)
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
    forml  <- paste( 'lmer(', target, '~', paste(pred_add, collapse = '+'), '+(1|site_a), data = NPP_statistical_final)')
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
  list_aic[[  a ]] <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = NPP_statistical_final)'))))
  
  list_bic[[ a ]] <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = NPP_statistical_final)'))))
  
  list_R[[ a ]] <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = NPP_statistical_final)'))))[1]
  preds_candidate <- preds_candidate[-which(preds_candidate == pred_max)]
}


R_null <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = NPP_statistical_final)'))))[1]
AIC_null <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = NPP_statistical_final)'))))
BIC_null <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = NPP_statistical_final)'))))
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

# the lowest AIC and BIC occurs at three factors (soil C/N + age + fAPAR)
p1;p2;p3
#####


#secondly
Nmin_statistical <- read.csv("/Users/yunpeng/data/NPP_final/Nmin_validation.csv")
Nmin_statistical <- subset(Nmin_statistical,Nmin>0)

Nmin_statistical$Nmin_a <- Nmin_statistical$Nmin

Nmin_statistical$age_a <- log(Nmin_statistical$age)
Nmin_statistical$alpha_a <- log(Nmin_statistical$alpha)
Nmin_statistical$Tg_a <- Nmin_statistical$Tg
Nmin_statistical$PPFD_a <- log(Nmin_statistical$PPFD)
Nmin_statistical$vpd_a <- log(Nmin_statistical$vpd)
Nmin_statistical$fAPAR_a <- Nmin_statistical$fAPAR
Nmin_statistical$CNrt_a <- log(Nmin_statistical$CNrt)
Nmin_statistical$LMA_a <- log(Nmin_statistical$LMA)
Nmin_statistical$vcmax25_a <- log(Nmin_statistical$max_vcmax25_c3)
Nmin_statistical$site_a <- Nmin_statistical$sitename

Nmin_statistical_large <- Nmin_statistical[,c("Nmin_a","age_a","alpha_a","Tg_a","PPFD_a","vpd_a","fAPAR_a",
                                            "CNrt_a","LMA_a","site_a")]

Nmin_statistical_final <- na.omit(Nmin_statistical_large)

#r2 = 0.66
summary(lmer(Nmin_a~Tg_a+alpha_a+(1|site_a),data=Nmin_statistical_final))


library(tidyverse)
library(ggplot2)

#1. select the most important predictor 
summary(Nmin_statistical_final)
#determine targets and preds.
target <- 'Nmin_a'

preds <- Nmin_statistical_final %>% select(-c(Nmin_a,site_a)) %>% 
  names()

r_list <- c()
library(lme4)
#For loop functions, include all predictor's r2 at the end
for (var in preds){
  forml <- paste( 'lmer(', target, '~', var, '+(1|site_a), data = Nmin_statistical_final)')
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

ggplot(All_rsquare, aes(x = reorder(preds, -rsq), y = rsq)) + geom_point() + theme(axis.text.x = element_text(angle = 60, hjust = 1))

#2. stepwise regression selection
library(caret)
library(recipes)
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
    forml  <- paste( 'lmer(', target, '~', paste(pred_add, collapse = '+'), '+(1|site_a), data = Nmin_statistical_final)')
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
  list_aic[[  a ]] <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = Nmin_statistical_final)'))))
  
  list_bic[[ a ]] <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = Nmin_statistical_final)'))))
  
  list_R[[ a ]] <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained, collapse = '+'),  '+(1|site_a), data = Nmin_statistical_final)'))))[1]
  preds_candidate <- preds_candidate[-which(preds_candidate == pred_max)]
}


R_null <- r.squaredGLMM(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = Nmin_statistical_final)'))))[1]
AIC_null <- AIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = Nmin_statistical_final)'))))
BIC_null <- BIC(eval(parse(text = paste( 'lmer(', target, '~', paste(preds_retained[1], collapse = '+'),  '+(1|site_a), data = Nmin_statistical_final)'))))
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

# the lowest AIC and BIC occurs at three factors (soil C/N + age + fAPAR)
p1;p2;p3


#now add additional map data
#check nc (unit: kg m-2 month-1) 1000*31556952/12
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")

#CABLE-POP
CABLE_POP_GPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CABLE-POP_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CABLE_POP_NPP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CABLE-POP_S2_npp_ANN_mean.nc"), varnam = "npp"))

#CABLE-POP
CABLE_POP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CABLE-POP_S2_gpp_ANN_mean.nc"), varnam = "gpp"))
CABLE_POP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CABLE-POP_S2_npp_ANN_mean.nc"), varnam = "npp"))


CABLE_POP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CABLE-POP_S2_npp_ANN_mean.nc"), varnam = "npp"))
CABLE_POP$myvar <- CABLE_POP$myvar#*1000*31556952/12

CLASS_CTEM <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CLASS-CTEM_S1_npp_ANN_mean.nc"), varnam = "npp"))
CLASS_CTEM$myvar <- CLASS_CTEM$myvar*1000*31556952/12

CLM5 <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/CLM5.0_S1_npp_ANN_mean.nc"), varnam = "npp"))
CLM5$myvar <- CLM5$myvar*1000*31556952/12
CLM5$lon[CLM5$lon>180] <- CLM5$lon[CLM5$lon>180]-360

ISAM <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ISAM_S1_npp_ANN_mean.nc"), varnam = "npp"))
ISAM$myvar <- ISAM$myvar*1000*31556952/12

JSBACH <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/JSBACH_S1_npp_ANN_mean.nc"), varnam = "npp"))
JSBACH$myvar <- JSBACH$myvar*1000*31556952/12

LPJ_GUESS <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/LPJ-GUESS_S1_npp_ANN_mean.nc"), varnam = "npp"))
LPJ_GUESS$myvar <- LPJ_GUESS$myvar*1000*31556952/12

ORCHIDEE <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ORCHIDEE_S1_npp_ANN_mean.nc"), varnam = "npp"))
ORCHIDEE$myvar <- ORCHIDEE$myvar*1000*31556952/12

#ORCHIDEE_CNP <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/ORCHIDEE-CNP_S1_npp_ANN_mean.nc"), varnam = "npp"))
#ORCHIDEE_CNP$myvar <- ORCHIDEE_CNP$myvar*1000*31556952/12

SDGVM <- as.data.frame(nc_to_df(read_nc_onefile("/Users/yunpeng/data/trendy/v8/SDGVM_S1_npp_ANN_mean.nc"), varnam = "npp"))
SDGVM$myvar <- SDGVM$myvar*1000*31556952/12


#
validation <- read.csv("~/data/NPP_final/NPP_validation.csv")
final_bnpp <- na.omit(validation[,c("lon","lat","TNPP_1")])

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


#resample --> small resolution
elev_nc <- read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc")
elev <- as.data.frame(nc_to_df(elev_nc, varnam = "elevation"))
names(elev) <- c("lon","lat","z")
coordinates(elev) <- ~lon+lat 
gridded(elev) <- TRUE
raster_half <- raster(elev, "z") 
bounding_box <- extent(-180, 180, -90, 90)
raster_half <- crop(raster_half, bounding_box)
raster_half

#resample --> big resolution --> only CABLE work
dim(na.omit(CABLE_POP))
dim(na.omit(CLASS_CTEM))
dim(na.omit(CLM5))
dim(na.omit(JSBACH))

CABLE_POP_df <- CABLE_POP
names(CABLE_POP_df) <- c("lon","lat","CABLE_POP")
coordinates(CABLE_POP_df) <- ~lon+lat 
gridded(CABLE_POP_df) <- TRUE
CABLE_POP_res <- raster(CABLE_POP_df, "CABLE_POP") 
CABLE_POP_res <- crop(CABLE_POP_res, bounding_box)

#others are error: > gridded(JSBACH_df) <- TRUE
#suggested tolerance minimum: 0.00844618 
#Error in points2grid(points, tolerance, round) : 
#dimension 2 : coordinate intervals are not constant

resampled_CABLE_POP <- raster::resample(CABLE_POP_res, raster_half, method="ngb")
CABLE_POP_final <- as.data.frame(stack(resampled_CABLE_POP),xy = TRUE)
names(CABLE_POP_final) <- c("lon","lat","CABLE_POP")


all_npp <- Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat"),all.x=TRUE),
                 list(na.omit(ORCHIDEE),na.omit(ISAM),na.omit(LPJ_GUESS),na.omit(SDGVM),na.omit(CABLE_POP_final),na.omit(all_maps[,c("lon","lat","npp_pft")])))

names(all_npp) <- c("lon","lat","ORCHIDEE","ISAM","LPJ_GUESS","SDGVM","CABLE_POP","Model")

#add alpha

alpha <- as.data.frame(nc_to_df(read_nc_onefile(
  "~/data/nimpl_sofun_inputs/map/Final_ncfile/alpha.nc"),
  varnam = "alpha"))
names(alpha) <- c("lon","lat","alpha")

final_npp <- Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat"),all.x=TRUE),
       list(all_npp,alpha,all_predictors))
       
model <- final_npp[,c("ORCHIDEE","ISAM","LPJ_GUESS","SDGVM","CABLE_POP","Model")]

predictors <- final_npp[,c("alpha","Tg","PPFD","vpd","fAPAR","age","CNrt","LMA","vcmax25")]

My_Theme = theme(
  axis.title.x = element_text(size = 15),
  axis.text.x = element_text(size = 20),
  axis.title.y = element_text(size = 15),
  axis.text.y = element_text(size = 20))

#alpha
for (a in 1:ncol(model)){
  gg <- heatscatter(predictors[,1],model[,a],xlab=colnames(predictors)[1],ylab=colnames(model)[a],ggplot=TRUE)+My_Theme+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=5)+My_Theme+xlim(0,1)#+geom_smooth(method = "lm", se = TRUE,color="blue")
  assign(paste0(colnames(predictors)[1],a), gg) 
}

plot_grid(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'),nrow=2,label_x = 0.9,label_y=0.92)
ggsave(paste("~/data/output/alpha1.jpg",sep=""),width = 20, height = 10)

#others
for (i in 2:ncol(predictors)){
  for (a in 1:ncol(model)){
    gg <- heatscatter(predictors[,i],model[,a],xlab=colnames(predictors)[i],ylab=colnames(model)[a],ggplot=TRUE)+My_Theme+
      stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 3,size=5)+My_Theme#+geom_smooth(method = "lm", se = TRUE,color="red")
    assign(paste0(colnames(predictors)[i],a), gg) 
  }
  print(i)
}

plot_grid(Tg1,Tg2,Tg3,Tg4,Tg5,Tg6,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'),nrow=2,label_x = 0.9,label_y=0.92)
ggsave(paste("~/data/output/Tg1.jpg",sep=""),width = 20, height = 10)

plot_grid(CNrt1,CNrt2,CNrt3,CNrt4,CNrt5,CNrt6,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'),nrow=2,label_x = 0.9,label_y=0.92)
ggsave(paste("~/data/output/CNrt1.jpg",sep=""),width = 20, height = 10)

plot_grid(age1,age2,age3,age4,age5,age6,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'),nrow=2,label_x = 0.9,label_y=0.92)
ggsave(paste("~/data/output/age1.jpg",sep=""),width = 20, height = 10)

#fapar
for (a in 1:ncol(model)){
  gg <- heatscatter(predictors[,5],model[,a],xlab=colnames(predictors)[5],ylab=colnames(model)[a],ggplot=TRUE)+My_Theme+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=5)+My_Theme+xlim(0,1)#+geom_smooth(method = "lm", se = TRUE,color="blue")
  assign(paste0(colnames(predictors)[5],a), gg) 
}

plot_grid(fAPAR1,fAPAR2,fAPAR3,fAPAR4,fAPAR5,fAPAR6,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'),nrow=2,label_x = 0.9,label_y=0.92)
ggsave(paste("~/data/output/fAPAR1.jpg",sep=""),width = 20, height = 10)


#finally, validation

validation <- read.csv("~/data/NPP_final/NPP_validation.csv")
final_bnpp <- na.omit(validation[,c("lon","lat","z","TNPP_1","pred_npp")])

plotinfo <- unique(final_bnpp[,c("lon","lat","z")])

a <- 1.5 # which degree (distance) of grid when interpolating gwr from global grids

model_xy <- final_npp[,c("lon","lat","ORCHIDEE","ISAM","LPJ_GUESS","SDGVM","CABLE_POP","Model")]
model_xyz <- Reduce(function(x,y) merge(x = x, y = y, by = c("lon","lat"),all.x=TRUE),
                    list(model_xy,all_maps[,c("lon","lat","z")]))

for (i in 1:nrow(plotinfo)) {
  tryCatch({
    #1
    Tg_global <- na.omit(model_xyz[,c("lon","lat","z","ORCHIDEE")])
    NRE_part <- subset(Tg_global,lon>(plotinfo[i,"lon"]-a)&lon<(plotinfo[i,"lon"]+a)&
                         lat>(plotinfo[i,"lat"]-a)&lat<(plotinfo[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- plotinfo[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    plotinfo[i,c("ORCHIDEE")] <- (gwr(ORCHIDEE ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #2
    Tg_global <- na.omit(model_xyz[,c("lon","lat","z","ISAM")])
    NRE_part <- subset(Tg_global,lon>(plotinfo[i,"lon"]-a)&lon<(plotinfo[i,"lon"]+a)&
                         lat>(plotinfo[i,"lat"]-a)&lat<(plotinfo[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- plotinfo[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    plotinfo[i,c("ISAM")] <- (gwr(ISAM ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #3
    Tg_global <- na.omit(model_xyz[,c("lon","lat","z","LPJ_GUESS")])
    NRE_part <- subset(Tg_global,lon>(plotinfo[i,"lon"]-a)&lon<(plotinfo[i,"lon"]+a)&
                         lat>(plotinfo[i,"lat"]-a)&lat<(plotinfo[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- plotinfo[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    plotinfo[i,c("LPJ_GUESS")] <- (gwr(LPJ_GUESS ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #4
    Tg_global <- na.omit(model_xyz[,c("lon","lat","z","SDGVM")])
    NRE_part <- subset(Tg_global,lon>(plotinfo[i,"lon"]-a)&lon<(plotinfo[i,"lon"]+a)&
                         lat>(plotinfo[i,"lat"]-a)&lat<(plotinfo[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- plotinfo[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    plotinfo[i,c("SDGVM")] <- (gwr(SDGVM ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
    #5
    Tg_global <- na.omit(model_xyz[,c("lon","lat","z","CABLE_POP")])
    NRE_part <- subset(Tg_global,lon>(plotinfo[i,"lon"]-a)&lon<(plotinfo[i,"lon"]+a)&
                         lat>(plotinfo[i,"lat"]-a)&lat<(plotinfo[i,"lat"]+a))
    coordinates(NRE_part) <- c("lon","lat")
    gridded(NRE_part) <- TRUE
    NRE_coord <- plotinfo[i,c("lon","lat","z")]
    coordinates(NRE_coord) <- c("lon","lat")
    plotinfo[i,c("CABLE_POP")] <- (gwr(CABLE_POP ~ z, NRE_part, bandwidth = 1.06, fit.points =NRE_coord,predictions=TRUE))$SDF$pred
  }, error=function(e){})} 

summary(plotinfo)
final_plot <- merge(final_bnpp,plotinfo,by=c("lon","lat","z"),all.x=TRUE)
names(final_plot)

a1 <- ggplot(data=final_plot, aes(x=pred_npp, y=TNPP_1)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  theme_classic()+My_Theme+labs(y = "Observed BP") +labs(x = "Our Predicted BP") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3)+xlim(0,2000)+ylim(0,2000)

a2 <- ggplot(data=final_plot, aes(x=ORCHIDEE, y=TNPP_1)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  theme_classic()+My_Theme+labs(y = "Observed BP") +labs(x = "ORCHIDEE BP") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3)+xlim(0,2000)+ylim(0,2000)

a3 <- ggplot(data=final_plot, aes(x=ISAM, y=TNPP_1)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  theme_classic()+My_Theme+labs(y = "Observed BP") +labs(x = "ISAM BP") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3)+xlim(0,2000)+ylim(0,2000)

a4 <- ggplot(data=final_plot, aes(x=LPJ_GUESS, y=TNPP_1)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  theme_classic()+My_Theme+labs(y = "Observed BP") +labs(x = "LPJ_GUESS BP") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3)+xlim(0,2000)+ylim(0,2000)

a5 <- ggplot(data=final_plot, aes(x=SDGVM, y=TNPP_1)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  theme_classic()+My_Theme+labs(y = "Observed BP") +labs(x = "SDGVM BP") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3)+xlim(0,2000)+ylim(0,2000)

a6 <- ggplot(data=final_plot, aes(x=CABLE_POP, y=TNPP_1)) +
  geom_point(alpha=0.5)+geom_abline(intercept=0,slope=1, linetype=3)+geom_smooth(method = "lm", se = F,size=2)+
  theme_classic()+My_Theme+labs(y = "Observed BP") +labs(x = "CABLE_POP BP") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3)+xlim(0,2000)+ylim(0,2000)

plot_grid(a2,a3,a4,a5,a6,a1,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'),nrow=2,label_x = 0.9,label_y=0.92)
ggsave(paste("~/data/output/bp_validation.jpg",sep=""),width = 20, height = 10)
