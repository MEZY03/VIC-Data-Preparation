##########################
# get_soil.py
##########################
import geopandas as gpd
import rasterio
import pandas as pd
import numpy as np
from rasterio.mask import mask
from rasterstats import zonal_stats
import matplotlib.pyplot as plt
from shapely.geometry import Point, box
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.family']=['Times New Roman','SimSun']
plt.rcParams['axes.unicode_minus'] = False

def extract_soil_to_grid(hwsd_raster, hwsd_csv, grid_coords_file, output_dir, grid_resolution=0.1):
	"""
	从HWSD中提取土壤数据到网格点 - 使用Excel导出的CSV文件
	
	参数:
	hwsd_raster: HWSD栅格文件路径 (hwsd.bil)
	hwsd_csv: 从Excel导出的HWSD_DATA.csv文件路径
	grid_coords_file: 网格坐标文件路径 (CSV格式)
	output_dir: 输出目录
	grid_resolution: 网格分辨率 (度), 如果为None则自动计算
	
	返回:
	soil_df: 土壤数据DataFrame
	"""
	
	print("=== HWSD土壤数据提取 (Excel CSV版本) ===")
	
	# 1. 检查输入文件
	print("1. 检查输入文件...")
	if not os.path.exists(hwsd_csv):
		print(f"错误: HWSD CSV文件不存在: {hwsd_csv}")
		print("请先用Excel导出HWSD_DATA表为CSV格式")
		return None
	if not os.path.exists(hwsd_raster):
		print(f"错误: HWSD栅格文件不存在: {hwsd_raster}")
		return None
	if not os.path.exists(grid_coords_file):
		print(f"错误: 网格坐标文件不存在: {grid_coords_file}")
		return None
	
	# 2. 读取网格坐标数据
	print("2. 读取网格坐标...")
	grid_coords = pd.read_csv(grid_coords_file)
	lons = grid_coords['lon_center']
	lats = grid_coords['lat_center']
	
	print(f"  网格点数量: {len(grid_coords)}")
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
	polygons = []
	for idx, row in grid_coords.iterrows():
		lon, lat = row['lon_center'], row['lat_center']
		half_res = grid_resolution / 2
		polygon = box(lon - half_res, lat - half_res, 
					 lon + half_res, lat + half_res)
		polygons.append(polygon)
	
	grid_gdf = gpd.GeoDataFrame(
		grid_coords,
		geometry=polygons,
		crs='EPSG:4326'
	)
	
	# 3. 读取HWSD CSV数据
	print("3. 读取HWSD CSV数据...")
	hwsd_data = pd.read_csv(hwsd_csv)
	print(f"HWSD数据记录数: {len(hwsd_data)}")
	print(f"HWSD数据列: {list(hwsd_data.columns)}")
	
	# 4. 计算每个网格的主导土壤类型
	print("4. 计算主导土壤类型...")
	with rasterio.open(hwsd_raster) as src:
		print(f"HWSD栅格范围: {src.bounds}")
		print(f"HWSD栅格尺寸: {src.width} x {src.height}")
		
		# 使用rasterstats计算众数
		stats = zonal_stats(
			grid_gdf, 
			hwsd_raster, 
			stats=['majority'], 
			categorical=True,
			nodata=src.nodata if src.nodata else -9999
		)
		
		majority_values = [stat.get('majority', np.nan) for stat in stats]
		grid_gdf['MU_GLOBAL'] = majority_values
		
		valid_count = sum(~pd.isna(majority_values))
		print(f"成功计算 {valid_count} 个网格的主导土壤类型")
		print(f"缺失土壤数据: {len(grid_gdf) - valid_count} 个网格")
	
	# 5. 连接土壤属性数据
	print("5. 连接土壤属性...")
	# 确保MU_GLOBAL类型一致
	grid_gdf['MU_GLOBAL'] = pd.to_numeric(grid_gdf['MU_GLOBAL'], errors='coerce')
	hwsd_data['MU_GLOBAL'] = pd.to_numeric(hwsd_data['MU_GLOBAL'], errors='coerce')
	
	# 执行连接
	merged_gdf = grid_gdf.merge(
		hwsd_data, 
		on='MU_GLOBAL', 
		how='left',
		suffixes=('', '_hwsd')
	)
	
	connection_rate = (sum(~merged_gdf['MU_GLOBAL'].isna()) / len(merged_gdf)) * 100
	print(f"土壤属性连接成功率: {connection_rate:.1f}%")
	
	# 6. 保存结果
	print("6. 保存结果...")
	os.makedirs(output_dir, exist_ok=True)
	
	# 保存土壤参数CSV文件
	output_csv = f"{output_dir}/grid_soil.csv"
	
	# 定义VIC模型需要的关键土壤参数
	vic_soil_parameters = [
		'grid_id', 'lon_center', 'lat_center', 'MU_GLOBAL',
		# 表层土壤参数 (Topsoil)
		'T_SAND', 'T_CLAY', 'T_SILT', 'T_OC',
		'T_GRAVEL', 'T_REF_BULK_DENSITY', 'T_ECE',
		# 底层土壤参数 (Subsoil)
		'S_SAND', 'S_CLAY', 'S_SILT', 'S_OC',
		'S_GRAVEL', 'S_REF_BULK_DENSITY', 'S_ECE'
	]
	
	# 只选择存在的列
	available_columns = [col for col in vic_soil_parameters if col in merged_gdf.columns]
	
	# 保存精简的CSV文件
	soil_df = merged_gdf[available_columns]
	soil_df.to_csv(output_csv, index=False, float_format='%.4f')
	print(f"土壤参数文件: {output_csv}")
	
	# 7. 生成可视化
	print("7. 生成可视化...")
	# 创建2行2列的子图布局
	fig, axes = plt.subplots(2, 2, figsize=(18, 12))
	
	# 子图1: 主导土壤类型分布
	valid_data = merged_gdf[~merged_gdf['MU_GLOBAL'].isna()]
	# 将土壤类型转换为整数
	valid_data['MU_GLOBAL_int'] = valid_data['MU_GLOBAL'].astype(int)
	valid_data.plot(column='MU_GLOBAL_int', ax=axes[0,0], legend=False,  # 如果土壤类型太多可关闭图例
					   cmap='tab20', categorical=True,
					   legend_kwds={'title': '土壤类型'})
	axes[0,0].set_title('主导土壤类型分布', fontsize=14)
	axes[0,0].set_xlabel('经度', fontsize=12)
	axes[0,0].set_ylabel('纬度', fontsize=12)
	axes[0,0].grid(True, alpha=0.3)
	# 子图位置微调
	current_pos = axes[0,0].get_position()
	axes[0,0].set_position([current_pos.x0-0.03, current_pos.y0, current_pos.width*0.97, current_pos.height*0.97])
	
	# 子图2: 表层砂粒含量分布
	valid_sand = merged_gdf[~merged_gdf['T_SAND'].isna()]
	valid_sand.plot(column='T_SAND', ax=axes[0,1], legend=True,
					   cmap='OrRd',
					   legend_kwds={'label': '砂粒含量 (%)'})
	axes[0,1].set_title('表层砂粒含量', fontsize=14)
	axes[0,1].set_xlabel('经度', fontsize=12)
	axes[0,1].set_ylabel('纬度', fontsize=12)
	axes[0,1].grid(True, alpha=0.3)
			
	# 子图3: 表层粘粒含量
	valid_clay = merged_gdf[~merged_gdf['T_CLAY'].isna()]
	valid_clay.plot(column='T_CLAY', ax=axes[1,0], legend=True,
					   cmap='YlGn',
					   legend_kwds={'label': '粘粒含量 (%)'})
	axes[1,0].set_title('表层粘粒含量', fontsize=14)
	axes[1,0].set_xlabel('经度', fontsize=12)
	axes[1,0].set_ylabel('纬度', fontsize=12)
	axes[1,0].grid(True, alpha=0.3)
	
	# 子图4: 表层有机碳含量
	valid_oc = merged_gdf[~merged_gdf['T_OC'].isna()]
	valid_oc.plot(column='T_OC', ax=axes[1,1], legend=True,
					   cmap='BuPu',
					   legend_kwds={'label': '有机碳含量 (%)'})
	axes[1,1].set_title('表层有机碳含量', fontsize=14)
	axes[1,1].set_xlabel('经度', fontsize=12)
	axes[1,1].set_ylabel('纬度', fontsize=12)
	axes[1,1].grid(True, alpha=0.3)
	
	vis_file = f"{output_dir}/soil_visualization.png"
	plt.savefig(vis_file, dpi=300, bbox_inches='tight')
	print(f"可视化文件: {vis_file}")
	
	plt.close()
	
	# 8. 输出详细的统计信息
	print("\n=== 土壤统计 ===")
	# 主导土壤类型统计
	if 'MU_GLOBAL' in merged_gdf.columns:
		soil_types = merged_gdf['MU_GLOBAL'].value_counts()
		print(f"  主导土壤类型数量: {len(soil_types)}")
		print(f"  前5种主要土壤类型:")
		for soil_type, count in soil_types.head().items():
			print(f"    MU_GLOBAL {int(soil_type)}: {count} 个网格 ({count/len(merged_gdf)*100:.1f}%)")
	
	# 关键土壤参数统计
	key_parameters = {
		'T_SAND': '表层砂粒含量',
		'T_CLAY': '表层粘粒含量', 
		'T_SILT': '表层粉粒含量',
		'T_OC': '表层有机碳含量',
		'T_GRAVEL': '表层砾石含量'
	}
	
	for param, desc in key_parameters.items():
		if param in merged_gdf.columns:
			valid_data = merged_gdf[param].dropna()
			if len(valid_data) > 0:
				print(f"  {desc}: {valid_data.min():.1f}-{valid_data.max():.1f}% "
					  f"(平均: {valid_data.mean():.1f}%)")
	
	print("\n=== 土壤数据提取完成 ===")
	return soil_df

def check_soil_texture(soil_df):
	"""
	检查VIC模型所需的土壤参数完整性
	"""
	print("\n=== VIC模型土壤参数检查 ===")
	
	# 检查VIC模型关键参数
	required_params = ['T_SAND', 'T_CLAY', 'T_SILT', 'T_OC']
	optional_params = ['T_GRAVEL', 'T_REF_BULK_DENSITY', 'S_SAND', 'S_CLAY', 'S_SILT']
	
	print("关键参数检查:")
	for param in required_params:
		if param in soil_df.columns:
			missing_rate = (soil_df[param].isna().sum() / len(soil_df)) * 100
			print(f"  {param}: {missing_rate:.1f}% 缺失")
		else:
			print(f"  {param}: 列不存在!")
	
	print("\n可选参数检查:")
	for param in optional_params:
		if param in soil_df.columns:
			missing_rate = (soil_df[param].isna().sum() / len(soil_df)) * 100
			print(f"  {param}: {missing_rate:.1f}% 缺失")
		else:
			print(f"  {param}: 列不存在")
	
	# 检查土壤质地总和
	if all(param in soil_df.columns for param in ['T_SAND', 'T_CLAY', 'T_SILT']):
		valid_rows = soil_df[['T_SAND', 'T_CLAY', 'T_SILT']].dropna()
		if len(valid_rows) > 0:
			total = valid_rows['T_SAND'] + valid_rows['T_CLAY'] + valid_rows['T_SILT']
			avg_total = total.mean()
			print(f"\n土壤质地验证:")
			print(f"  砂粒+粘粒+粉粒平均总和: {avg_total:.1f}%")
			if abs(avg_total - 100) > 5:
				print(f"  警告: 土壤质地总和偏离100%较多，可能需要检查数据")

# 完整使用流程示例
if __name__ == "__main__":
	# 设置文件路径
	hwsd_raster = "hwsd.bil"  # HWSD栅格文件
	hwsd_csv = "HWSD_DATA.csv"  # 从Excel导出的CSV文件
	grid_coords_file = "../grid/vic_grid/grid_coordinates.csv"  # 网格坐标CSV
	output_dir = "./soil_data"  # 输出目录
	grid_resolution = 0.1  # VIC网格分辨率
	
	# 提取土壤数据
	soil_df = extract_soil_to_grid(
		hwsd_raster=hwsd_raster,
		hwsd_csv=hwsd_csv,
		grid_coords_file=grid_coords_file,
		output_dir=output_dir,
		grid_resolution=grid_resolution
	)
	
	# 输出结果摘要
	print(f"\n输出文件:")
	print(f"  - 土壤数据: {output_dir}/grid_soil.csv")
	print(f"  - 可视化图表: {output_dir}/soil_visualization.png")
	
	# 检查VIC模型参数完整性
	check_soil_texture(soil_df)
	
	# 显示前几个网格的土壤参数
	print(f"\n前5个网格的土壤参数:")
	print(soil_df.head())