-----Preprocessing------

---Global prediction fields
1. c3c4: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/c3c4.Rmd
2. PPFD + vpd + Tg + alpha: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/climates.Rmd
3. fAPAR + age: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/fAPAR_age.Rmd
4. LMA: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/LMA.Rmd
5. Soil C/N: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/soilCN.Rmd
6. vcmax25: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/Vcmax25.Rmd


---Prepare complete dataset for measurements:
7. forest npp + grassland npp + Nuptake: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/Forest_site_orig.R


#please note! 8-10 below requires csv that output from 7 ("~/data/NPP_Yunke/NPP_Nmin_dataset.csv"). This csv in the future may be changed/updated. However, as far as lon, lat, z, begin_year and end_year not changed/added. It's fine to not replicate the (too long) process below.

---Forcing and p-model: (but not used in New Phytol submission since we don't simulate GPP here)
8. Simulated GPP for forest npp + grassland npp + Nuptake: pmodel_simulation.R. This takes time, and output data in ~/data/NPP_Yunke/simulated_gpp/site_simulated_gpp_vcmax.csv

---collect site level Tg, vpd, alpha and PPFD basing on measurement year, using gwr. 
9. This take time and see climate_site_data.R. This takes time, and output data in ...

---collect site level Tg, vpd, alpha and PPFD from Map, using gwr. 
10. This take time and see map_site_data.R.


---Now, output figures in manuscript!

12. Statistical model: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/statistical_fig1_tab1.R
save all statsitical model info in: file = "~/data/NPP_Grassland_final/statistical_model/...RData"
	
13. valdation figure: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/validation_fig2.R

14. Global simulation: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/global_figs.R

15: Uncertainty: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission/uncertainty_table_s1.R

16: AIC for model selection: 
