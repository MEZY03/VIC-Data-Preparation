# 修改脚本权限
# find ./ -type d -exec chmod 755 {} \;
# find ./ -type f -exec chmod 644 {} \;
find ./ -name "*.sh" -exec chmod 777 {} \;

# 运行创建网格脚本
cd grid
python create_grid.py
cd ..

# 运行创建高程脚本
cd elevation
python get_elevation.py
cd ..

# 运行创建降水脚本
cd precipitation
./calculate_climatology.sh
python get_precipitation.py
cd ..

# 运行创建土壤参数脚本
cd soil
python get_soil.py
python calculate_soil_water.py
python create_soil_param.py
cd ..

# 运行创建植被参数脚本
cd veg
python process_modis.py
python create_veg_param.py
cd ..

# 运行创建强迫数据脚本
cd forcing
rm -r forcings
./get_forcing.sh
python create_forcing.py
cd ..

# 整合并打包所有结果
cp -f soil/soil_param.txt veg/veg_param.txt veg/veg_lib_IGBP.txt vic_input_parameters
tar -czvf forcings.tar.gz -C forcing forcings
(cd vic_input_parameters && tar -czvf ../parameters.tar.gz *.txt)