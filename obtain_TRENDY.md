#below is an exmple for obtain TRENDY
ssh yunpeng@euler.ethz.ch

cd /cluster/work/climate/bestocke/data/trendy/v8

sftp trendy-v8@trendy.ex.ac.uk

gcp-2019

cd output

cd CABLE-POP

cd S2

get CABLE-POP_S2_npp.nc

cd ../../

exit 

module load gcc/4.8.2 cdo/1.6.4

#all others
cdo selyear,1986/2014 LPX-Bern_S2_gpp.nc a1.nc
cdo muldpm a1.nc a2.nc
cdo mulc,86400000 a2.nc a3.nc
cdo yearsum a3.nc a4.nc
cdo -O selyear,1986/2015 a4.nc a5.nc
cdo -O timmean a5.nc LPX-Bern_S2_gpp_ANN_mean.nc

#DLEM --> cannot conver to normal nc

#ISAM
cdo seltimestep,127/156 ISAM_S2_fNup.nc a1.nc
#cdo mulc,1000 -vertsum a1.nc a2.nc
#cdo seltimestep,127/156 a2.nc a3.nc
cdo -O timmean a4.nc ISAM_S2_fNup_ANN_mean.nc


CLM5.0_S2_fNup.nc
