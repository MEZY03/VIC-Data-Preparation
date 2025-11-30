##########################
# create_grid.py
##########################
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import box
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.family']=['Times New Roman','SimSun']
plt.rcParams['axes.unicode_minus'] = False

def create_vic_grid(boundary_file, output_dir, grid_resolution):
	"""
	为VIC模型创建规则网格系统
	
	参数:
	boundary_file: 研究区域边界文件路径 (Shapefile)
	output_dir: 输出目录路径
	grid_resolution: 网格分辨率 (度)
	
	返回:
	grid_data: 包含网格几何和中心坐标的数据GeoDataFrame
	"""
	
	print("=== VIC模型网格生成 ===")
	
	# 1. 读取边界文件
	print("1. 读取研究区域边界...")
	boundary = gpd.read_file(boundary_file)
	
	# 转换为WGS84坐标系
	if boundary.crs != 'EPSG:4326':
		boundary = boundary.to_crs('EPSG:4326')
		print("  已转换为WGS84坐标系")
	
	# 计算网格参数
	bounds = boundary.total_bounds
	min_lon, min_lat, max_lon, max_lat = bounds
	
	print(f"  研究区域范围: {min_lon:.4f}, {min_lat:.4f} - {max_lon:.4f}, {max_lat:.4f}")
	print(f"  网格分辨率: {grid_resolution}度")
	
	n_cols = int((max_lon - min_lon) / grid_resolution) + 1
	n_rows = int((max_lat - min_lat) / grid_resolution) + 1
	print(f"  网格尺寸: {n_cols}列 × {n_rows}行")
	
	# 2. 生成完整网格
	print("2. 生成网格系统...")
	grid_polygons = []
	grid_centers = []
	grid_ids = []
	
	grid_id = 0
	for i in range(n_rows):
		for j in range(n_cols):
			# 计算网格边界
			left = min_lon + j * grid_resolution
			right = left + grid_resolution
			bottom = min_lat + i * grid_resolution
			top = bottom + grid_resolution
			
			# 创建网格多边形和中心点
			polygon = box(left, bottom, right, top)
			center_lon = (left + right) / 2
			center_lat = (bottom + top) / 2
			
			grid_polygons.append(polygon)
			grid_centers.append((center_lon, center_lat))
			grid_ids.append(grid_id)
			grid_id += 1
	
	# 创建完整网格数据框
	full_grid = gpd.GeoDataFrame({
		'grid_id': grid_ids,
		'geometry': grid_polygons,
		'lon_center': [p[0] for p in grid_centers],
		'lat_center': [p[1] for p in grid_centers]
	}, crs='EPSG:4326')
	
	# 3. 裁剪到研究区域
	print("3. 裁剪网格到研究区域...")
	clipped_grid = gpd.sjoin(full_grid, boundary, how='inner', predicate='intersects')
	grid_data = clipped_grid[['grid_id', 'geometry', 'lon_center', 'lat_center']]
	
	print(f"  完整网格数: {len(full_grid)}")
	print(f"  研究区域内网格数: {len(clipped_grid)}")
	
	# 4. 检查并处理重复网格
	print("4. 检查数据完整性...")
	duplicate_count = grid_data.duplicated(subset=['grid_id']).sum()
	
	if duplicate_count > 0:
		print(f"  发现 {duplicate_count} 个重复网格ID，进行去重处理")
		grid_data = grid_data.drop_duplicates(subset=['grid_id'], keep='first')
		print(f"  去重后网格数: {len(grid_data)}")
	else:
		print("  网格ID唯一性检查通过")
	
	# 5. 保存结果文件
	print("5. 保存输出文件...")
	# 创建输出目录
	os.makedirs(output_dir, exist_ok=True)
	
	# 保存网格坐标CSV
	grid_csv = f"{output_dir}/grid_coordinates.csv"
	grid_data[['grid_id', 'lon_center', 'lat_center']].to_csv(grid_csv, index=False)
	print(f"  网格坐标文件: {grid_csv}")
	
	# 保存VIC格式网格信息
	vic_info = f"{output_dir}/vic_grid_info.txt"
	with open(vic_info, 'w') as f:
		f.write("# VIC模型网格信息\n")
		f.write("# grid_id lat lon\n")
		for _, row in grid_data.iterrows():
			f.write(f"{row['grid_id']} {row['lat_center']:.6f} {row['lon_center']:.6f}\n")
	print(f"  VIC网格信息: {vic_info}")
	
	# 6. 生成可视化
	print("6. 生成可视化图表...")
	fig, ax = plt.subplots(figsize=(12, 8))
	
	# 绘制研究区域边界
	boundary.boundary.plot(ax=ax, color='red', linewidth=2, label='研究区域')
	
	# 绘制网格
	grid_data.boundary.plot(ax=ax, color='blue', linewidth=0.5, alpha=0.7, label='VIC网格')
	
	# 绘制网格中心点
	ax.scatter(grid_data['lon_center'], grid_data['lat_center'], 
			   color='green', s=10, alpha=0.5, label='网格中心')
	
	ax.set_title(f'VIC模型网格系统 ({len(grid_data)}个网格)', fontsize=14)
	ax.set_xlabel('经度', fontsize=12)
	ax.set_ylabel('纬度', fontsize=12)
	ax.legend()
	ax.grid(True, alpha=0.3)
	
	# 保存图片
	vis_file = f"{output_dir}/grid_visualization.png"
	plt.tight_layout()
	plt.savefig(vis_file, dpi=300, bbox_inches='tight')
	print(f"  可视化文件: {vis_file}")
	
	plt.close()
	
	print("\n=== 网格生成完成 ===")
	print(f"成功创建 {len(grid_data)} 个网格单元")
	
	return grid_data

# 完整使用流程示例
if __name__ == "__main__":
	# 参数设置
	boundary_file = "Lower_YRB.shp"      # 研究区域边界文件
	output_dir = "./vic_grid"      # 输出目录
	grid_resolution = 0.1                # 网格分辨率 (度)
	
	# 执行网格生成
	grid_data = create_vic_grid(
		boundary_file=boundary_file,
		output_dir=output_dir,
		grid_resolution=grid_resolution
	)
	
	# 输出结果摘要
	print(f"\n输出文件:")
	print(f"  - 网格坐标文件: {output_dir}/grid_coordinates.csv")
	print(f"  - VIC网格信息: {output_dir}/vic_grid_info.txt")
	print(f"  - 可视化图表: {output_dir}/grid_visualization.png")
	
	# 显示前几个网格信息
	print(f"\n前5个网格信息:")
	print(grid_data[['grid_id', 'lon_center', 'lat_center']].head())