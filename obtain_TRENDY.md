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

#all others (except for LPX-Bern that uses 1986-2014, since in obtain_TRENDY_year.sh it claims that 2015 has some issue, all others using 1986-2015)
cdo selyear,1901/2015 CLM5.0_S2_fNup.nc a1.nc # a1 is monthly data for kgC/m2/s
cdo muldpm a1.nc a2.nc # a2 is monthly data for kgC/m2/s * 30 or 31 (converted from day to month)
cdo mulc,86400000 a2.nc a3.nc # a3 monthly data for convert from seconds to day, and unit from kg to g --> now a3 is gC/m2/month
cdo yearsum a3.nc a4.nc  # a4 is yearly data convert from gC/m2/month to gC/m2/year
cdo -O selyear,1986/2015 a4.nc a5.nc #select year again - select years again so that the value become stable.
cdo -O timmean a5.nc CLM5.0_S2_fNup_ANN_mean.nc # calculate average of these years

#DLEM --> cannot conver to normal nc

#ISAM --> unit still in kg/m2/s --> needs to convert to g/m2/yr in R code (by multiplying 1000*31556952) then
cdo seltimestep,127/156 ISAM_S2_fNup.nc a1.nc
cdo -O timmean a1.nc ISAM_S2_fNup_ANN_mean.nc