# CDO预处理代码
cdo -houravg -select,name=TLML     MERRA/M2T1NXFLX.5.12.4_MERRA2_*.tavg1_2d_flx_Nx.*.nc4.nc4 forcing_data/temperature.nc
cdo -houravg -select,name=PS       MERRA/M2T1NXSLV.5.12.4_MERRA2_*.tavg1_2d_slv_Nx.*.nc4.nc4 forcing_data/pressure.nc
cdo -houravg -select,name=QLML     MERRA/M2T1NXFLX.5.12.4_MERRA2_*.tavg1_2d_flx_Nx.*.nc4.nc4 forcing_data/specific_humidity.nc
cdo -houravg -select,name=SPEED    MERRA/M2T1NXFLX.5.12.4_MERRA2_*.tavg1_2d_flx_Nx.*.nc4.nc4 forcing_data/wind_speed.nc
cdo -houravg -select,name=SWGDN    MERRA/M2T1NXRAD.5.12.4_MERRA2_*.tavg1_2d_rad_Nx.*.nc4.nc4 forcing_data/shortwave_radiation.nc
cdo -houravg -select,name=LWGAB    MERRA/M2T1NXRAD.5.12.4_MERRA2_*.tavg1_2d_rad_Nx.*.nc4.nc4 forcing_data/longwave_radiation.nc
cdo -houravg -select,name=PRECTOT  MERRA/M2T1NXFLX.5.12.4_MERRA2_*.tavg1_2d_flx_Nx.*.nc4.nc4 forcing_data/precipitation.nc

# 如果变量维顺序错误, 需要用以下NCO代码调换
# ncpdq --rdr=time,lat,lon input_file.nc output_file.nc