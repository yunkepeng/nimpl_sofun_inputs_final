-----Preprocessing

---Global prediction fields
1. c3c4: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/c3c4.Rmd
2. PPFD + vpd + Tg + alpha: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/climates.Rmd
3. fAPAR + age: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/fAPAR_age.Rmd
4. LMA: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/LMA.Rmd
5. Soil C/N: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/soilCN.Rmd
6. vcmax25: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/Vcmax25.Rmd


---Prepare Validation data for forest and grassland...

7. write final dataset of forest + grassland csv (including finally checked sitename for forcing and rep info):
/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_site_orig.R

output two csv in ("/Users/yunpeng/data/NPP_final/NPP_Forest.csv") and ("/Users/yunpeng/data/NPP_Grassland_final/NPP_grassland.csv") and 

---Prepare Training data for forest

8. prepare training data for forest statistical model: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/NPP_statistical_model.Rmd
output: 
csvfile <- paste("/Users/yunpeng/data/NPP_final/NPP_statistical_forest.csv")
csvfile <- paste("/Users/yunpeng/data/NPP_final/NRE_statistical_forest.csv")


---Forcing and p-model 

9. Forest simulation: nimpl_sofun_inputs_final/NPP/Forest_simulation.R 

Already finsished for calculating all predicted c and n uptakes, but I have repeated again in validation file
output:
csvfile <- paste("/Users/yunpeng/data/NPP_final/NPP_validation.csv")
csvfile <- paste("/Users/yunpeng/data/NPP_final/Nmass_validation.csv")
csvfile <- paste("/Users/yunpeng/data/NPP_final/NRE_validation.csv")

Also output siteinfo in "/Users/yunpeng/data/NPP_final/fpar_name/forest_fpar_name.csv" (just note to myself!)

10. Grassland simulation: nimpl_sofun_inputs_final/NPP/Grassland_simulation.R
FINISH by line 577!!!
Already finsished for calculating all predicted c and n uptakes, but I have repeated again in validation file
output:
csvfile <- paste("/Users/yunpeng/data/NPP_Grassland_final/NPP_grass_validation.csv")

Also output siteinfo in "/Users/yunpeng/data/NPP_final/fpar_name/grassland_fpar_name.csv" (just note to myself!)


11. N minerlization simulation: nimpl_sofun_inputs_final/NPP/New_Nuptake_site_simulation.R
FINISH by line 426!!!
output:
csvfile <- paste("/Users/yunpeng/data/NPP_final/Nmin_validation.csv")

---Now, output figures in manuscript!

12. Statistical model: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/statistical_fig1_tab1.R
save all statsitical model info in: file = "~/data/NPP_Grassland_final/statistical_model/...RData"
	
13. valdation figure: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/validation_fig2.R

14. Global simulation: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/global_figs.R

15: Uncertainty: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/uncertainty_table_s1.R

16: AIC for model selection: 
