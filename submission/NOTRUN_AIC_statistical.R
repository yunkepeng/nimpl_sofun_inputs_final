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
#check AIC of each models, using stepwise regression design

###Data preparation 
NPP_statistical <- read.csv("~/data/NPP_final/NPP_statistical_forest.csv")
NPP_statistical$tnpp_gpp_a <- log((NPP_statistical$TNPP_1/NPP_statistical$GPP)/(1-(NPP_statistical$TNPP_1/NPP_statistical$GPP)))
NPP_statistical$soilCN_a <- log(NPP_statistical$soilCN)
NPP_statistical$age_a <- log(NPP_statistical$age)
NPP_statistical$observedfAPAR_a <- NPP_statistical$observedfAPAR
NPP_statistical$alpha_a <- log(NPP_statistical$alpha_gwr)
NPP_statistical$PPFD_a <- log(NPP_statistical$PPFD_gwr)
NPP_statistical$Tg_a <- NPP_statistical$Tg_gwr
NPP_statistical$vpd_a <- log(NPP_statistical$vpd_gwr)
NPP_statistical$site_a <- NPP_statistical$site
NPP_statistical_final <- NPP_statistical[,c("tnpp_gpp_a","soilCN_a","age_a","observedfAPAR_a",
                                            "alpha_a","PPFD_a","Tg_a","vpd_a","site_a")]
#################################
#1. Now, prepare Stepwise regression for NPP/GPP
#################################
library(tidyverse)
library(ggplot2)

#1. select the most important predictor 
summary(NPP_statistical_final)
#determine targets and preds.
target <- 'tnpp_gpp_a'

preds <- NPP_statistical_final %>% select(-c(tnpp_gpp_a,site_a)) %>% 
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


#################################
#2. Now, prepare Stepwise regression for ANPP/GPP
#################################
NPP_statistical$anpp_gpp_a <- log((NPP_statistical$ANPP_2/NPP_statistical$GPP)/(1-(NPP_statistical$ANPP_2/NPP_statistical$GPP)))
NPP_statistical_final <- NPP_statistical[,c("anpp_gpp_a","soilCN_a","age_a","observedfAPAR_a",
                                            "alpha_a","PPFD_a","Tg_a","vpd_a","site_a")]

#1. select the most important predictor 
summary(NPP_statistical_final)
#determine targets and preds.
target <- 'anpp_gpp_a'

preds <- NPP_statistical_final %>% select(-c(anpp_gpp_a,site_a)) %>% 
  names()

r_list <- c()
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
p4 <- ggplot() + 
  geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = R_all)) 
#AIC
p5 <- ggplot() + 
  geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = AIC_all)) 
#BIC
p6 <- ggplot() + 
  geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = BIC_all)) 

#the lowest AIC and BIC occurs at 4 factors (soil C/N + age + fAPAR + PPFD), while additoinal factor of PPFD only increases from R2 from 0.41 to 0.43 (not necessary), so removed.
p4;p5;p6


#################################
#3. Now, prepare Stepwise regression for leaf-NPP/ANPP
#################################
#first - inluding all 6 factors like above two models (soil C/N + age+ fAPAR + climates) - but the best model was only when including soil C/N and the sign is unexpected positive - so we alternatively to use a climate-driven model
NPP_statistical$lnpp_anpp_a <- log((NPP_statistical$NPP.foliage/NPP_statistical$ANPP_2)/(1-(NPP_statistical$NPP.foliage/NPP_statistical$ANPP_2)))
NPP_statistical$lnpp_anpp_a[NPP_statistical$NPP.foliage>=NPP_statistical$ANPP_2] <- NA
lnpp_data <- subset(NPP_statistical,file=="Sara Vicca"|file=="ForC")
NPP_statistical_final <-NPP_statistical[,c("lnpp_anpp_a","soilCN_a","age_a","observedfAPAR_a",
                                                                     "PPFD_a","Tg_a","vpd_a","site_a")]
NPP_statistical_final <- na.omit(NPP_statistical_final)

#determine targets and preds.
target <- 'lnpp_anpp_a'

preds <- NPP_statistical_final %>% select(-c(lnpp_anpp_a,site_a)) %>% 
  names()

r_list <- c()
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
p7 <- ggplot() + 
  geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = R_all)) 
#AIC
p8 <- ggplot() + 
  geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = AIC_all)) 
#BIC
p9 <- ggplot() + 
  geom_point(data = df1, aes(x = factor(variable_all,level = variable_all), y = BIC_all)) 

summary(lmer(lnpp_anpp_a ~ age_a  + (1|site_a),data=NPP_statistical_final))
#the lowst AIC and BIC occured when only including soil C/N and the sign is unexpected positive - so we alternatively to use a climate-driven model
p7;p8;p9
summary(lmer(lnpp_anpp_a~Tg_a+vpd_a+PPFD_a+(1|site_a),data=lnpp_data))
r.squaredGLMM(lmer(lnpp_anpp_a~Tg_a+vpd_a+PPFD_a+(1|site_a),data=lnpp_data))
