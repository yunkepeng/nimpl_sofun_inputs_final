#Forest NPP validation: calculations see ~/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
NPP_validation <- read.csv("~/data/NPP_final/NPP_validation.csv")
NPP_validation$obs_age[NPP_validation$obs_age==999] <- NA

summary(NPP_validation$TNPP_1)
#measured variable: soilCN,observedfAPAR,obs_age
#Peredicted variable: age, alpha, Tg,PPFD,vpd, fAPAR, CNrt,LMA,max_vcmax25_c3


NPP_statistical <- read.csv("~/data/NPP_final/NPP_validation.csv")

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



#check nc
elev_nc <- read_nc_onefile("/Users/yunpeng/data/trendy/output_mean.nc")
#elev_nc <- read_nc_onefile("D:/PhD/nimpl_sofun_inputs/Data/Elevation/WFDEI-elevation.nc")
elev <- as.data.frame(nc_to_df(elev_nc, varnam = "npp"))
summary(elev*1000000000)
