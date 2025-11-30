##########################
# get_precipitation.py
##########################
import geopandas as gpd
import xarray as xr
import numpy as np
import pandas as pd
from shapely.geometry import Point, box
import matplotlib.pyplot as plt
import os
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.family']=['Times New Roman','SimSun']
plt.rcParams['axes.unicode_minus'] = False

def extract_precipitation_to_grid(precip_file, grid_coords_file, output_dir, precip_variable_name, grid_resolution=None):
	"""
	从降水气候态文件中提取降水数据到网格点
	
	参数:
	precip_file: 降水气候态NetCDF文件路径
	grid_coords_file: 网格坐标文件路径 (CSV格式)
	output_dir: 输出目录路径
	precip_variable_name: NetCDF中的降水变量名
	grid_resolution: 网格分辨率 (度), 如果为None则自动计算
	
	返回:
	precip_df: 降水数据DataFrame
	"""
	
	print("=== 提取降水数据到网格点 ===")
	
	# 1. 检查输入文件
	print("1. 检查输入文件...")
	if not os.path.exists(precip_file):
		print(f"错误: 降水气候态NetCDF文件不存在: {precip_file}")
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
	
	# 3. 读取降水气候态数据
	print("3. 读取降水气候态数据...")
	try:
		ds = xr.open_dataset(precip_file)
		
		# 检查变量名
		if precip_variable_name not in ds.variables:
			available_vars = list(ds.variables.keys())
			print(f"变量 '{precip_variable_name}' 不存在")
			print(f"  可用变量: {available_vars}")
			# 尝试使用第一个数据变量
			data_vars = [var for var in available_vars if len(ds[var].shape) >= 2]
			if data_vars:
				precip_variable_name = data_vars[0]
				print(f"  使用变量: {precip_variable_name}")
			else:
				raise KeyError("未找到合适的数据变量")
		
		precip_data = ds[precip_variable_name]
		
		# 获取坐标信息
		coord_mapping = {
			'longitude': ['longitude', 'lon', 'x'],
			'latitude': ['latitude', 'lat', 'y']
		}
		
		lons, lats = None, None
		for std_name, possible_names in coord_mapping.items():
			for name in possible_names:
				if name in ds.coords:
					if std_name == 'longitude':
						lons = ds[name].values
					else:
						lats = ds[name].values
					break
		
		if lons is None or lats is None:
			# 从数据维度推断
			dims = precip_data.dims
			if len(dims) >= 2:
				lons = ds[dims[-1]].values if dims[-1] in ds.coords else None
				lats = ds[dims[-2]].values if dims[-2] in ds.coords else None
		
		if lons is None or lats is None:
			raise ValueError("无法确定经纬度坐标")
		
		print(f"  数据维度: {precip_data.dims}")
		print(f"  数据形状: {precip_data.shape}")
		print(f"  经度范围: {lons.min():.2f} ~ {lons.max():.2f}")
		print(f"  纬度范围: {lats.min():.2f} ~ {lats.max():.2f}")
		print(f"  降水范围: {precip_data.min().values:.1f} ~ {precip_data.max().values:.1f} mm/year")
		print(f"  数据单位: {precip_data.attrs.get('units', '未知')}")
		
		# 获取数据值
		precip_values = precip_data.values
		
		# 处理二维网格数据
		if len(precip_values.shape) == 2:
			lon_grid, lat_grid = np.meshgrid(lons, lats)
			points_2d_source = np.column_stack([lon_grid.ravel(), lat_grid.ravel()])
			values_source = precip_values.ravel()
		else:
			raise ValueError(f"不支持的降水数据维度: {precip_values.shape}")
			
		ds.close()
		
	except Exception as e:
		print(f"读取降水数据失败: {e}")
		return None
	
	# 4. 插值气候态降水到网格点
	print("4. 插值到网格点...")
	
	target_points = np.column_stack([grid_gdf['lon_center'], grid_gdf['lat_center']])
	
	precip_interpolated = griddata(
		points_2d_source, 
		values_source, 
		target_points, 
		method='linear',
		fill_value=np.nan
	)
	
	grid_gdf['precipitation'] = precip_interpolated
	
	valid_count = np.sum(~np.isnan(precip_interpolated))
	missing_count = len(grid_gdf) - valid_count
	print(f"  成功插值 {valid_count} 个网格点")
	print(f"  缺失气候态降水数据: {missing_count} 个点")
	if missing_count > 0:
		print("  注意: 部分网格点位于降水数据范围之外")
	
	# 5. 保存结果
	print("5. 保存结果...")
	os.makedirs(output_dir, exist_ok=True)
	
	# 保存降水数据CSV文件
	output_csv = f"{output_dir}/grid_precipitation.csv"
	precip_df = grid_gdf[['grid_id', 'lon_center', 'lat_center', 'precipitation']]
	precip_df.to_csv(output_csv, index=False)
	print(f"  降水数据文件: {output_csv}")
	
	# 6. 生成可视化
	print("6. 生成可视化图表...")
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
	
	# 子图1: 网格降水分布
	valid_data = grid_gdf[~grid_gdf['precipitation'].isna()]
	valid_data.plot(column='precipitation', ax=ax1, legend=True,
					   cmap='Blues',marker='s',
					   legend_kwds={'label': '年降水量 (mm/year)'})

	ax1.set_title('网格降水分布', fontsize=14)
	ax1.set_xlabel('经度', fontsize=12)
	ax1.set_ylabel('纬度', fontsize=12)
	ax1.grid(True, alpha=0.3)
	
	# 子图2: 降水统计直方图
	valid_precip = grid_gdf['precipitation'].dropna()
	ax2.hist(valid_precip, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
	ax2.axvline(valid_precip.mean(), color='red', linestyle='--', linewidth=2, 
			   label=f'平均值: {valid_precip.mean():.1f} mm/year')
	ax2.set_title('降水值分布', fontsize=14)
	ax2.set_xlabel('年降水量 (mm/year)', fontsize=12)
	ax2.set_ylabel('频数', fontsize=12)
	ax2.legend()
	ax2.grid(True, alpha=0.3)
	
	plt.tight_layout()
	vis_file = f"{output_dir}/precipitation_visualization.png"
	plt.savefig(vis_file, dpi=300, bbox_inches='tight')
	print(f"  可视化文件: {vis_file}")
	
	plt.close()
	
	# 7. 输出详细的统计信息
	print("\n=== 降水统计 ===")
	if len(valid_precip) > 0:
		print(f"  最小值: {valid_precip.min():.1f} mm/year")
		print(f"  最大值: {valid_precip.max():.1f} mm/year")
		print(f"  平均值: {valid_precip.mean():.1f} mm/year")
		print(f"  标准差: {valid_precip.std():.1f} mm/year")
		print(f"  中位数: {valid_precip.median():.1f} mm/year")
		print(f"  缺失率: {missing_count/len(grid_gdf)*100:.1f}%")
		
		# 空间分布统计
		print(f"  空间分布:")
		print(f"    - 最高降水区域: {grid_gdf.loc[grid_gdf['precipitation'].idxmax(), ['lon_center', 'lat_center']].values}")
		print(f"    - 最低降水区域: {grid_gdf.loc[grid_gdf['precipitation'].idxmin(), ['lon_center', 'lat_center']].values}")
	else:
		print("  无有效降水数据")
	
	print("\n=== 降水数据提取完成 ===")
	return precip_df

# 完整使用流程示例
if __name__ == "__main__":
	# 参数设置
	precip_file = "./precipitation_data/precip_climatology.nc"  # 降水气候态文件
	grid_coords_file = "../grid/vic_grid/grid_coordinates.csv"  # 网格坐标文件
	output_dir = "./precipitation_data"                   # 输出目录
	precip_variable_name = "PRECTOT"  # 降水变量名
	grid_resolution = 0.1  # VIC网格分辨率
	
	# 执行降水提取
	precip_df = extract_precipitation_to_grid(
		precip_file=precip_file,
		grid_coords_file=grid_coords_file,
		output_dir=output_dir,
		precip_variable_name=precip_variable_name,
		grid_resolution=grid_resolution
	)
	
	# 输出结果摘要
	print(f"\n输出文件:")
	print(f"  - 降水数据: {output_dir}/grid_precipitation.csv")
	print(f"  - 可视化图表: {output_dir}/precipitation_visualization.png")
	
	# 显示前几个网格的降水信息
	print(f"\n前5个网格降水:")
	print(precip_df.head())