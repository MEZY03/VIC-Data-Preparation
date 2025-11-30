##########################
# create_veg_param.py
##########################
import geopandas as gpd
import pandas as pd
import numpy as np
from osgeo import gdal
from rasterstats import zonal_stats
import os
import warnings
warnings.filterwarnings('ignore')

def calculate_vegetation_fractions(soil_param_file, landcover_file, grid_resolution=None):
	"""
	计算各植被类型的比例
	
	参数:
	soil_param_file: 土壤参数文件路径
	landcover_file: 土地覆盖栅格文件路径
	grid_resolution: 网格分辨率 (度), 如果为None则自动计算
	
	返回:
	vegetation_data: 植被比例数据list
	"""
	
	print("=== 开始计算植被类型比例 ===")
	
	# 1. 从土壤参数文件读取网格信息
	print("1. 读取网格信息...")
	soil_params = pd.read_csv(soil_param_file, delim_whitespace=True, header=None)
	
	# 提取网格ID、纬度和经度
	grid_ids = soil_params[1].values  # 第2列是grid_id
	lats = soil_params[2].values      # 第3列是纬度
	lons = soil_params[3].values      # 第4列是经度
	
	print(f"  网格点数量: {len(grid_ids)}")
	print(f"  经度范围: {lons.min():.3f} ~ {lons.max():.3f}")
	print(f"  纬度范围: {lats.min():.3f} ~ {lats.max():.3f}")
	
	# 处理网格大小
	if grid_resolution is None:
		if len(lons) > 1:
			# 使用经度差自动计算网格大小 (如果未提供)
			grid_resolution = abs(lons[1] - lons[0])
			print(f"  自动计算的空间分辨率: {grid_resolution:.6f} 度")
		else:
			# 如果只有一个网格，使用默认值
			grid_resolution = 0.1
			print(f"  警告: 只有一个网格，使用默认分辨率: {grid_resolution} 度")
	else:
		print(f"  使用指定的空间分辨率: {grid_resolution:.3f} 度")
	
	# 创建网格多边形
	grid_features = []
	
	for i, (grid_id, lat, lon) in enumerate(zip(grid_ids, lats, lons)):
		half_res = grid_resolution / 2
		left = lon - half_res
		right = lon + half_res
		bottom = lat - half_res
		top = lat + half_res
		
		polygon = {
			'type': 'Polygon',
			'coordinates': [[
				[left, bottom],
				[right, bottom], 
				[right, top],
				[left, top],
				[left, bottom]
			]]
		}
		
		grid_features.append({
			'type': 'Feature',
			'properties': {'grid_id': grid_id, 'lat': lat, 'lon': lon},
			'geometry': polygon
		})
	
	# 转换为GeoDataFrame
	grid_gdf = gpd.GeoDataFrame.from_features(grid_features)
	grid_gdf.crs = "EPSG:4326"
	
	# 2. IGBP到VIC植被类型映射
	print("2. 映射IGBP到VIC...")
	IGBP_TO_VIC = {
		1: 1,    # Evergreen Needleleaf Forest
		2: 2,    # Evergreen Broadleaf Forest
		3: 3,    # Deciduous Needleleaf Forest
		4: 4,    # Deciduous Broadleaf Forest
		5: 5,    # Mixed Forests
		6: 6,    # Closed Shrublands
		7: 7,    # Open Shrublands
		8: 8,    # Woody Savannas
		9: 9,    # Savannas
		10: 10,  # Grasslands
		11: 11,  # Permanent Wetlands
		12: 12,  # Croplands
		13: 13,  # Urban and Built-up
		14: 14,  # Cropland/Natural Vegetation Mosaic
		15: 15,  # Snow and Ice
		16: 16,  # Barren
		17: 0,   # Water (IGBP 17)
	}
	
	# 3. 执行区域统计
	print("3. 正在进行区域统计...")
	stats = zonal_stats(grid_gdf, 
					   landcover_file,
					   categorical=True,
					   category_map=IGBP_TO_VIC,
					   all_touched=True,
					   nodata=255)
	
	# 4. 处理统计结果
	print("4. 处理统计结果...")
	vegetation_data = []
	for i, (grid_id, lat, lon, stat) in enumerate(zip(grid_ids, lats, lons, stats)):
		total_pixels = sum(stat.values())
		
		if total_pixels == 0:
			print(f"警告: 网格 {grid_id} 没有有效数据")
			# 添加默认植被类型（裸地）
			vegetation_data.append({
				'grid_id': grid_id,
				'lat': lat,
				'lon': lon,
				'vegetation_fractions': [(16, 1.0)]  # 100% 裸地
			})
			continue
		
		# 计算植被比例 (各像素数除以总像素数)
		vic_fractions = {}
		for veg_type, pixel_count in stat.items():
			vic_contribution = pixel_count / total_pixels
			if veg_type in vic_fractions:
				vic_fractions[veg_type] += vic_contribution
			else:
				vic_fractions[veg_type] = vic_contribution
		
		# 归一化处理
		total_fraction = sum(vic_fractions.values())
		if abs(total_fraction - 1.0) > 0.01:
			for vic_type in vic_fractions:
				vic_fractions[vic_type] /= total_fraction
		
		sorted_fractions = sorted(vic_fractions.items(), key=lambda x: x[1], reverse=True)
		
		vegetation_data.append({
			'grid_id': grid_id,
			'lat': lat, 
			'lon': lon,
			'vegetation_fractions': sorted_fractions
		})
	
	print(f"植被比例计算完成，处理了 {len(vegetation_data)} 个网格")
	return vegetation_data


def generate_veg_param(soil_param_file, vegetation_data, output_file, veg_lai_params=None, veg_albedo_params=None):
	"""
	基于 Schaperow (2021) 表格生成VIC模型所需的植被参数文件
	
	参数:
	soil_param_file: 土壤参数文件路径
	vegetation_data: 之前计算的植被比例数据
	output_file: 输出文件路径
	veg_lai_params: LAI参数，如果为None则不输出LAI
	veg_albedo_params: 反照率参数，如果为None则不输出反照率
	"""
	
	print("=== 开始生成植被参数文件 ===")
	
	# 1. 检查参数设置
	if veg_lai_params is None:
		print("注意: veg_lai_params=None, 将不在veg_param.txt中输出LAI")
		print("提示: 请在全局参数文件中设置 LAI_SRC=FROM_VEGLIB")
	if veg_albedo_params is None:
		print("注意: veg_albedo_params=None, 将不在veg_param.txt中输出反照率")
		print("提示: 请在全局参数文件中设置 ALB_SRC=FROM_VEGLIB")
	
	# 2. 从土壤参数文件读取网格信息
	soil_params = pd.read_csv(soil_param_file, delim_whitespace=True, header=None)
	
	# 提取网格ID、纬度和经度
	grid_ids = soil_params[1].values  # 第2列是grid_id
	lats = soil_params[2].values      # 第3列是纬度
	lons = soil_params[3].values      # 第4列是经度
	
	print(f"从土壤参数文件中处理 {len(grid_ids)} 个网格单元")
	
	# 3. 植被类型参数默认值
	# 根区参数 - 假设3个根区
	veg_root_params = {
		# veg_type: [root_depth1, root_fract1, root_depth2, root_fract2, root_depth3, root_fract3]
		0:  [0.1, 0.44, 0.6, 0.45, 0.8, 0.11],   # 水体 (Open water)
		1:  [0.1, 0.34, 0.6, 0.51, 1.1, 0.14],   # 常绿针叶林 (Evergreen needleleaf forest)
		2:  [0.1, 0.32, 0.6, 0.44, 2.3, 0.23],   # 常绿阔叶林 (Evergreen broadleaf forest)
		3:  [0.1, 0.34, 0.6, 0.50, 1.3, 0.16],   # 落叶针叶林 (Deciduous needleleaf forest)
		4:  [0.1, 0.31, 0.6, 0.52, 1.3, 0.17],   # 落叶阔叶林 (Deciduous broadleaf forest)
		5:  [0.1, 0.25, 0.6, 0.52, 1.7, 0.22],   # 混交林 (Mixed forest)
		6:  [0.1, 0.31, 0.6, 0.49, 1.8, 0.21],   # 郁闭灌丛 (Closed shrublands)
		7:  [0.1, 0.33, 0.6, 0.43, 2.4, 0.24],   # 开放灌丛 (Open shrublands)
		8:  [0.1, 0.37, 0.6, 0.50, 1.0, 0.13],   # 木本稀树草原 (Woody savanna)
		9:  [0.1, 0.36, 0.6, 0.45, 1.7, 0.19],   # 稀树草原 (Savanna)
		10: [0.1, 0.44, 0.6, 0.45, 0.8, 0.11],   # 草地 (Grasslands)
		11: [0.1, 0.44, 0.6, 0.45, 0.8, 0.11],   # 湿地 (Permanent wetlands)
		12: [0.1, 0.33, 0.6, 0.55, 0.8, 0.12],   # 耕地 (Cropland)
		13: [0.1, 0.44, 0.6, 0.45, 0.8, 0.11],   # 城市 (Urban)
		14: [0.1, 0.33, 0.6, 0.55, 0.8, 0.12],   # 耕地/自然混合 (Cropland/natural vegetation mosaic)
		15: [0.1, 0.44, 0.6, 0.45, 0.8, 0.11],   # 冰雪 (Permanent snow and ice)
		16: [0.1, 0.22, 0.6, 0.46, 3.3, 0.31]    # 裸地 (Barren)
	}
	
	# 4. LAI月值 (根据表格数据更新)
	if veg_lai_params is None:
		# 使用空字典，表示不输出LAI
		veg_lai_params = {}
	elif veg_lai_params=='default':
		# 使用提供的LAI参数
		veg_lai_params = {
			# veg_type: 12个月的叶面积指数 (1月到12月)
			0:  [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],  # 水体 (Open Water)
			1:  [3.40, 3.40, 3.50, 3.70, 4.00, 4.40, 4.40, 4.30, 4.20, 3.70, 3.50, 3.40],  # 常绿针叶林 (Evergreen Needleleaf forest)
			2:  [3.40, 3.40, 3.50, 3.70, 4.00, 4.40, 4.40, 4.30, 4.20, 3.70, 3.50, 3.40],  # 常绿阔叶林 (Evergreen Broadleaf forest)
			3:  [1.68, 1.52, 1.68, 2.90, 4.90, 5.00, 5.00, 4.60, 3.44, 3.04, 2.16, 2.00],  # 落叶针叶林 (Deciduous Needleleaf forest)
			4:  [1.68, 1.52, 1.68, 2.90, 4.90, 5.00, 5.00, 4.60, 3.44, 3.04, 2.16, 2.00],  # 落叶阔叶林 (Deciduous Broadleaf forest)
			5:  [1.68, 1.52, 1.68, 2.90, 4.90, 5.00, 5.00, 4.60, 3.44, 3.04, 2.16, 2.00],  # 混交林 (Mixed forest)
			6:  [2.00, 2.25, 2.95, 3.85, 3.75, 3.50, 3.55, 3.20, 3.30, 2.85, 2.60, 2.20],  # 郁闭灌丛 (Closed Shrublands)
			7:  [2.00, 2.25, 2.95, 3.85, 3.75, 3.50, 3.55, 3.20, 3.30, 2.85, 2.60, 2.20],  # 开放灌丛 (Open Shrublands)
			8:  [1.68, 1.52, 1.68, 2.90, 4.90, 5.00, 5.00, 4.60, 3.44, 3.04, 2.16, 2.00],  # 木本稀树草原 (Woody Savannas)
			9:  [2.00, 2.25, 2.95, 3.85, 3.75, 3.50, 3.55, 3.20, 3.30, 2.85, 2.60, 2.20],  # 稀树草原 (Savannas)
			10: [2.00, 2.25, 2.95, 3.85, 3.75, 3.50, 3.55, 3.20, 3.30, 2.85, 2.60, 2.20],  # 草地 (Grasslands)
			11: [2.00, 2.25, 2.95, 3.85, 3.75, 3.50, 3.55, 3.20, 3.30, 2.85, 2.60, 2.20],  # 湿地 (Permanent wetlands)
			12: [0.50, 0.50, 0.50, 0.50, 1.50, 3.00, 4.50, 5.00, 2.50, 0.50, 0.50, 0.50],  # 耕地 (Croplands)
			13: [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],  # 城市 (Urban and built-up)
			14: [2.00, 2.25, 2.95, 3.85, 3.75, 3.50, 3.55, 3.20, 3.30, 2.85, 2.60, 2.20],  # 耕地/自然混合 (Cropland/Natural vegetation mosaic)
			15: [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],  # 冰雪 (Snow and ice)
			16: [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20]   # 裸地 (Barren or sparsely vegetated)
		}

	# 5. 反照率月值 (根据表格数据更新)
	if veg_albedo_params is None:
		# 使用空字典，表示不输出反照率
		veg_albedo_params = {}
	elif veg_albedo_params=='default':
		# 使用提供的反照率参数
		veg_albedo_params = {
			# veg_type: 12个月的反照率 (1月到12月)
			0:  [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],  # 水体 (Open Water)
			1:  [0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12],  # 常绿针叶林 (Evergreen Needleleaf forest)
			2:  [0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12],  # 常绿阔叶林 (Evergreen Broadleaf forest)
			3:  [0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18],  # 落叶针叶林 (Deciduous Needleleaf forest)
			4:  [0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18],  # 落叶阔叶林 (Deciduous Broadleaf forest)
			5:  [0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18],  # 混交林 (Mixed forest)
			6:  [0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19],  # 郁闭灌丛 (Closed Shrublands)
			7:  [0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19],  # 开放灌丛 (Open Shrublands)
			8:  [0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18],  # 木本稀树草原 (Woody Savannas)
			9:  [0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19],  # 稀树草原 (Savannas)
			10: [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],  # 草地 (Grasslands)
			11: [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],  # 湿地 (Permanent wetlands)
			12: [0.10, 0.10, 0.10, 0.10, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10],  # 耕地 (Croplands)
			13: [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20],  # 城市 (Urban and built-up)
			14: [0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19],  # 耕地/自然混合 (Cropland/Natural vegetation mosaic)
			15: [0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60],  # 冰雪 (Snow and ice)
			16: [0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20]   # 裸地 (Barren or sparsely vegetated)
		}
	
	# 6. 生成植被参数文件
	with open(output_file, 'w') as f:
		for i, (grid_id, lat, lon) in enumerate(zip(grid_ids, lats, lons)):
			# 获取该网格的植被类型比例数据
			if i < len(vegetation_data):
				grid_veg_data = vegetation_data[i]
			else:
				print(f"警告: 网格 {grid_id} 没有植被数据，跳过")
				continue
			
			# 网格号和植被类型数量
			n_veg_types = len(grid_veg_data['vegetation_fractions'])
			f.write(f"{int(grid_id)} {n_veg_types}\n")
			
			# 遍历该网格的每种植被类型
			for veg_type, fraction in grid_veg_data['vegetation_fractions']:
				veg_type = int(veg_type)
				
				# 植被类型和覆盖比例
				f.write(f"  {veg_type} {fraction:.4f}   ")
				
				# 根区参数（3个根区）
				if veg_type in veg_root_params:
					root_params = veg_root_params[veg_type]
					for zone in range(3):
						if zone==2:
							f.write(f"{root_params[zone*2]:.3f} {root_params[zone*2+1]:.3f}\n")
						else:
							f.write(f"{root_params[zone*2]:.3f} {root_params[zone*2+1]:.3f} ")
				else:
					# 使用默认值
					for zone in range(3):
						if zone==2:
							f.write("0.100 0.333\n")
						else:
							f.write("0.100 0.333 ")
				
				# LAI月值（12个月）- 仅在veg_lai_params不为None时输出
				if veg_lai_params and veg_type in veg_lai_params:
					lai_values = veg_lai_params[veg_type]
					lai_line = " ".join([f"{val:.3f}" for val in lai_values])
					f.write(f"    {lai_line}\n")
				elif veg_lai_params:
					# 使用默认值
					lai_line = " ".join(["1.000"] * 12) + "\n"
					f.write(f"    {lai_line}\n")
				# 如果veg_lai_params为None，则不输出LAI行
				
				# 反照率月值（12个月）- 仅在veg_albedo_params不为None时输出
				if veg_albedo_params and veg_type in veg_albedo_params:
					albedo_values = veg_albedo_params[veg_type]
					albedo_line = " ".join([f"{val:.3f}" for val in albedo_values])
					f.write(f"    {albedo_line}\n")
				elif veg_albedo_params:
					# 使用默认值
					albedo_line = " ".join(["0.150"] * 12) + "\n"
					f.write(f"    {albedo_line}\n")
				# 如果veg_albedo_params为None，则不输出反照率行
	
	print(f"VIC植被参数文件已生成: {output_file}")


# 完整使用流程示例
if __name__ == "__main__":
	# 参数设置
	soil_param_file = "../soil/soil_param.txt"  # 土壤参数文件
	landcover_file = "./modis_data/MCD12Q1_LC_Type1_mosaic.tif"  # 土地覆盖数据
	output_veg_file = "veg_param.txt"  # 输出植被参数文件
	grid_resolution = 0.1  # VIC网格分辨率
	
	# 计算植被比例
	vegetation_data = calculate_vegetation_fractions(soil_param_file, landcover_file, grid_resolution)
	
	# 生成植被参数文件
	generate_veg_param(soil_param_file, vegetation_data, output_veg_file, veg_lai_params='default', veg_albedo_params='default')
	
	# 输出结果摘要
	print(f"\n输出文件:")
	print(f"  - 植被参数: {output_veg_file}")
	
	# 显示基本网格统计信息
	n_grids = len(vegetation_data)
	print(f"网格总数: {n_grids}")
	
	# 每个网格的植被类型数量统计
	veg_types_per_grid = [len(grid['vegetation_fractions']) for grid in vegetation_data]
	avg_veg_types = np.mean(veg_types_per_grid)
	max_veg_types = np.max(veg_types_per_grid)
	min_veg_types = np.min(veg_types_per_grid)
	print(f"每个网格的植被类型数:")
	print(f"  - 平均: {avg_veg_types:.2f}")
	print(f"  - 最大: {max_veg_types}")
	print(f"  - 最小: {min_veg_types}")