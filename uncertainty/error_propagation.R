#error propagation 
rm(list=ls())
load("/Users/yunpeng/data/gpp_gmd/new/out_eval_FULL.Rdata")
obs <-out_eval_FULL$gpp$fluxnet$data$meandf$obs
pred <- out_eval_FULL$gpp$fluxnet$data$meandf$mod
plot(obs~pred)
obs_pred <- as.data.frame(cbind(obs,pred))
obs_pred <- na.omit(obs_pred)
summary(obs_pred)
# According to textbook:Uncertainty estimates obtained as standard deviations of repeated measurement results are called A type uncertainty estimates. 
# In this way,  “sample mean” corresponds to the true GPP, and the “observation” is the modelled GPP.

obs_pred$variance <- (obs_pred$obs - obs_pred$pred)^2
uncertainty_gpp <- sqrt(sum(obs_pred$variance)/nrow(obs_pred))
uncertainty_gpp


