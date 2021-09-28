---Global prediction fields
1. c3c4: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/Prediction_field/c3c4.Rmd
2. PPFD + vpd + Tg + alpha (not used): /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/Prediction_field/climates.Rmd
3. fAPAR + age: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/Prediction_field/fAPAR_age.Rmd
4. LMA: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/Prediction_field/LMA.Rmd
5. Soil C/N: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/Prediction_field/soilCN.Rmd
6. vcmax25: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/Prediction_field/Vcmax25.Rmd

All saved in: ~/data/nimpl_sofun_inputs/map/Final_ncfile/eg.nc

---Prepare Validation data for forest and grassland
7. write final dataset of forest + grassland csv (including finally checked sitename for forcing and rep info):
/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_site_orig.R

output two csv in ("/Users/yunpeng/data/NPP_final/NPP_Forest.csv") and ("/Users/yunpeng/data/NPP_Grassland_final/NPP_grassland.csv") and 

---Prepare Training data for forest
8. prepare training data for forest statistical model: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/NPP_statistical_model.Rmd
output: 
csvfile <- paste("/Users/yunpeng/data/NPP_final/NPP_statistical_forest.csv")
csvfile <- paste("/Users/yunpeng/data/NPP_final/NRE_statistical_forest.csv")


---Forcing and p-model 
9. Forest simulation: nimpl_sofun_inputs_final/NPP/Forest_simulation.Rmd 

Already finsished for calculating all predicted c and n uptakes, but I have repeated again in validation file
output:
csvfile <- paste("/Users/yunpeng/data/NPP_final/NPP_validation.csv")
csvfile <- paste("/Users/yunpeng/data/NPP_final/Nmass_validation.csv")
csvfile <- paste("/Users/yunpeng/data/NPP_final/NRE_validation.csv")

10. Grassland simulation: nimpl_sofun_inputs_final/NPP/Grassland_simulation.Rmd 
FINISH by line 577!!!
Already finsished for calculating all predicted c and n uptakes, but I have repeated again in validation file
output:
csvfile <- paste("/Users/yunpeng/data/NPP_Grassland_final/NPP_grass_validation.csv")

---Now, output figures in manuscript!
11. Statistical model: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/statistical_fig1_tab1.R
save all statsitical model info in: file = "~/data/NPP_Grassland_final/statistical_model/...RData"

12. valdation figure: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/validation_fig2.R

13. Global simulation: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/global_figs.R

14: Uncertainty: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/uncertainty_table_s1.R

15: AIC for model selection: 
