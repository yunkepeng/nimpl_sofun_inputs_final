Constant were saved in 

###constant of leaf C (%) in forest = 46% (mean of combined_leaf_traits.csv)
/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/NPP_statistical_model.Rmd
L206-208

SP_input <- read.csv(file="/Users/yunpeng/data/leaf_traits/combined_leaf_traits_updated.csv") #new one
summary(SP_input$C_percent)
#mean of value of leaf C = 0.46 or 46%



###constant of wood and root cn in forest (median):
/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/Forest_simulation.R
L543-549

#median of wood cn and root cn
# see below
summary(read.csv("/Users/yunpeng/data/CN_wood/wood_cn.csv")$OrigValueStr) #from TRY database
summary(NPP_Forest2$CN_root_final)
#using median of wood =100
#using median of root = 94




###constant of leaf and root cn in grassland (median):
/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/Grassland_simulation.R
L510-514

#leaf c/n model. median = 18
summary(NPP_grassland_final5_gpp_npp_anpp$CN_leaf_final)
#root c/n model median = 41
summary(NPP_grassland_final5_gpp_npp_anpp$CN_root_final)



###constant of NRE in grassland
/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/NPP/NPP_statistical_model.Rmd
L486-487 (last line)
summary(NRE_df$NRE)
#median of NRE in grassland is 0.69
