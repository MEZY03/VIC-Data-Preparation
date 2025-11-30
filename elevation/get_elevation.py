##########################
# get_elevation.py
##########################
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import Point, box
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.family']=['Times New Roman','SimSun']
plt.rcParams['axes.unicode_minus'] = False

def extract_elevation_to_grid(dem_file, grid_coords_file, output_dir, grid_resolution=None):
	"""
	从DEM数据提取网格点高程
	
	参数:
	dem_file: DEM栅格文件路径
	grid_coords_file: 网格坐标文件路径 (CSV格式)
	output_dir: 输出目录路径
	grid_resolution: 网格分辨率 (度), 如果为None则自动计算
	
	返回:
	elevation_df: 高程数据DataFrame
	"""
	
	print("=== 高程数据提取 ===")
	
	# 1. 检查输入文件
	print("1. 检查输入文件...")
	if not os.path.exists(dem_file):
		print(f"错误: DEM栅格文件不存在: {dem_file}")
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
	
	# 3. 读取DEM数据
	print("3. 读取DEM数据...")
	with rasterio.open(dem_file) as src:
		dem_crs = src.crs
		dem_transform = src.transform
		dem_data = src.read(1)
		dem_nodata = src.nodata if src.nodata is not None else -9999
		
		print(f"  DEM范围: {src.bounds}")
		print(f"  DEM分辨率: {src.res}")
		print(f"  DEM坐标系: {src.crs}")
	
	# 4. 提取网格点高程
	print("4. 提取高程值...") 
	# 创建点几何对象并转换坐标系
	points = [Point(lon, lat) for lon, lat in zip(grid_gdf['lon_center'], grid_gdf['lat_center'])]
	points_gdf = gpd.GeoDataFrame(grid_gdf, geometry=points, crs='EPSG:4326')
	points_gdf_proj = points_gdf.to_crs(dem_crs)
	
	elevations = []
	valid_count = 0
	for idx, row in points_gdf_proj.iterrows():
		x, y = row.geometry.x, row.geometry.y
		
		# 地理坐标转像素坐标
		col, row_idx = ~dem_transform * (x, y)
		col, row_idx = int(col), int(row_idx)
		
		# 检查坐标有效性并提取高程
		if 0 <= row_idx < dem_data.shape[0] and 0 <= col < dem_data.shape[1]:
			elevation = dem_data[row_idx, col]
			if elevation != dem_nodata and not np.isnan(elevation):
				elevations.append(elevation)
				valid_count += 1
			else:
				elevations.append(np.nan)
		else:
			elevations.append(np.nan)
	
	grid_gdf['elevation'] = elevations
	
	missing_count = len(grid_gdf) - valid_count
	print(f"  成功提取 {valid_count} 个点的高程数据")
	print(f"  缺失高程数据: {missing_count} 个点")
	
	# 5. 保存结果
	print("5. 保存结果...")
	os.makedirs(output_dir, exist_ok=True)
	
	# 保存高程数据CSV文件
	output_csv = f"{output_dir}/grid_elevation.csv"
	elevation_df = grid_gdf[['grid_id', 'lon_center', 'lat_center', 'elevation']]
	elevation_df.to_csv(output_csv, index=False)
	print(f"  高程数据文件: {output_csv}")
	
	# 6. 生成可视化
	print("6. 生成可视化图表...")
	fig, ax = plt.subplots(figsize=(10, 8))

	# 使用geopandas的plot方法
	valid_data = grid_gdf[~grid_gdf['elevation'].isna()]
	valid_data.plot(column='elevation', ax=ax, legend=True,
					   cmap='terrain',marker='s',
					   legend_kwds={'label': '高程 (m)'})
	ax.set_title('网格高程分布', fontsize=14)
	ax.set_xlabel('经度', fontsize=12)
	ax.set_ylabel('纬度', fontsize=12)
	ax.grid(True, alpha=0.3)

	plt.tight_layout()
	vis_file = f"{output_dir}/elevation_visualization.png"
	plt.savefig(vis_file, dpi=300, bbox_inches='tight')
	print(f"  可视化文件: {vis_file}")
	
	plt.close()
	
	# 7. 输出详细的统计信息
	print("\n=== 高程统计 ===")
	valid_elevations = grid_gdf['elevation'].dropna()
	if len(valid_elevations) > 0:
		print(f"  最小值: {valid_elevations.min():.1f} m")
		print(f"  最大值: {valid_elevations.max():.1f} m")
		print(f"  平均值: {valid_elevations.mean():.1f} m")
		print(f"  标准差: {valid_elevations.std():.1f} m")
	
	print("\n=== 高程数据提取完成 ===")
	return elevation_df

# 完整使用流程示例
if __name__ == "__main__":
	# 参数设置
	dem_file = "output_hh.tif"                # DEM文件路径
	grid_coords_file = "../grid/vic_grid/grid_coordinates.csv"  # 网格坐标文件
	output_dir = "./elevation_data"     # 输出目录
	grid_resolution = 0.1  # VIC网格分辨率
	
	# 执行高程提取
	elevation_df = extract_elevation_to_grid(
		dem_file=dem_file,
		grid_coords_file=grid_coords_file,
		output_dir=output_dir,
		grid_resolution=grid_resolution
	)
	
	# 输出结果摘要
	print(f"\n 输出文件:")
	print(f"  - 高程数据: {output_dir}/grid_elevation.csv")
	print(f"  - 可视化图表: {output_dir}/elevation_visualization.png")
	
	# 显示前几个网格的高程信息
	print(f"\n 前5个网格高程:")
	print(elevation_df.head())