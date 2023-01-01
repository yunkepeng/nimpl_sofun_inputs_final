-----Preprocessing------

---Global prediction fields
1. c3c4: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/c3c4.Rmd
2. PPFD + vpd + Tg + alpha: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/climates.Rmd
3. fAPAR + age: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/fAPAR_age.Rmd
4. LMA: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/LMA.Rmd
5. Soil C/N: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/soilCN.Rmd
6. vcmax25: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/Vcmax25.Rmd


---Prepare complete dataset for measurements:
7. forest npp + grassland npp + Nuptake + NRE: /Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/preprocessing/Forest_site_orig.R

#final dataset of NPP, Nuptake saved in /Users/yunpeng/data/NPP_Yunke/NPP_Nmin_dataset_with_predictors.csv
#final dataset of NRE saved in /Users/yunpeng/data/NRE_various/NRE_dataset.csv
#final data of leaf traits saved in /Users/yunpeng/data/leaf_traits/combined_leaf_traits_updated.csv (from Peng et al. 2021 Communications Biology)

#please note! 8-10 below requires csv that output from 7 ("~/data/NPP_Yunke/NPP_Nmin_dataset.csv"). This csv in the future may be changed/updated. However, as far as lon, lat, z, begin_year and end_year not changed/added. It's fine to not replicate the (too long) process below.

---Forcing and p-model: (but not used in New Phytol submission since we don't simulate GPP here)
8. Simulated GPP for forest npp + grassland npp + Nuptake: pmodel_simulation.R. This takes time, and output data in ~/data/NPP_Yunke/simulated_gpp/site_simulated_gpp_vcmax.csv

---collect site level Tg, vpd, alpha and PPFD basing on measurement year, using gwr. 
9. This take time and see climate_site_data.R. This takes time, and output data in ~data/NPP_Yunke/predictors/climates_sites.csv

---collect site level Tg, vpd, alpha and PPFD from Map, using gwr. 
10. This take time and see map_site_data.R.


---Now,my analysis
/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/submission_newphytol/submission/all_analyses.R