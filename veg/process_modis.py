##########################
# process_modis.py
##########################
from osgeo import gdal
import numpy as np
import os
from glob import glob
import warnings
warnings.filterwarnings('ignore')

def extract_modis_hdf(hdf_files, output_dir):
	"""
	从MCD12Q1 HDF文件中提取土地覆盖数据层
	
	参数:
	hdf_files: HDF文件列表
	output_dir: 输出目录路径
	
	返回:
	extracted_files: 提取的TIFF文件路径列表list
	"""
	layers = {
		'LC_Type1': 'HDF4_EOS:EOS_GRID:"{fname}":MCD12Q1:LC_Type1',
		'QC': 'HDF4_EOS:EOS_GRID:"{fname}":MCD12Q1:QC'
	}
	
	extracted_files = []
	
	for hdf_file in hdf_files:
		base_name = os.path.basename(hdf_file).replace('.hdf', '')
		
		for layer_name, template in layers.items():
			gdal_path = template.format(fname=hdf_file)
			output_tif = os.path.join(output_dir, f"{base_name}_{layer_name}.tif")
			
			try:
				gdal.Translate(output_tif, gdal_path)
				extracted_files.append(output_tif)
				print(f"提取完成: {output_tif}")
			except Exception as e:
				print(f"提取失败 {output_tif}: {e}")
	
	return extracted_files

def reproject_and_mosaic(input_files, output_dir, target_epsg="EPSG:4326", output_res=0.005):
	"""
	将正弦投影转换为地理坐标系并拼接
	
	参数:
	input_files: 输入TIFF文件列表
	output_dir: 输出目录路径
	target_epsg: 目标坐标系, 默认EPSG:4326
	output_res: 输出分辨率, 默认0.005度
	
	返回:
	mosaic_output: 拼接后的文件路径str
	"""
	reprojected_files = []
	
	for tif_file in input_files:
		if "LC_Type1" in tif_file:
			base_name = os.path.basename(tif_file).replace('.tif', '')
			output_file = os.path.join(output_dir, f"{base_name}_reprojected.tif")
			
			gdal.Warp(output_file, tif_file, 
					 dstSRS=target_epsg,
					 resampleAlg=gdal.GRA_NearestNeighbour,
					 xRes=output_res, yRes=output_res)
			reprojected_files.append(output_file)
			print(f"重投影完成: {output_file}")
	
	# 拼接
	mosaic_output = os.path.join(output_dir, "MCD12Q1_LC_Type1_mosaic.tif")
	
	if len(reprojected_files) > 1:
		gdal.Warp(mosaic_output, reprojected_files, 
				 format='GTiff', resampleAlg=gdal.GRA_NearestNeighbour)
		print(f"拼接完成: {mosaic_output}")
	else:
		gdal.Translate(mosaic_output, reprojected_files[0])
		print(f"单文件转换完成: {mosaic_output}")
	
	return mosaic_output

def clip_to_shp_region(input_raster, region_shp, output_dir):
	"""
	将数据裁剪到研究区域范围
	
	参数:
	input_raster: 输入栅格文件路径
	region_shp: 区域边界Shapefile路径
	output_dir: 输出目录路径
	
	返回:
	clipped_output: 裁剪后的文件路径str
	"""
	base_name = os.path.basename(input_raster).replace('.tif', '')
	clipped_output = os.path.join(output_dir, f"{base_name}_Yangtze.tif")
	
	gdal.Warp(clipped_output, input_raster,
			  cutlineDSName=region_shp,
			  cropToCutline=True,
			  dstNodata=255,
			  resampleAlg=gdal.GRA_NearestNeighbour)
	
	print(f"裁剪完成: {clipped_output}")
	return clipped_output

def process_modis_data(hdf_files, output_dir, region_shp=None):
	"""
	处理MODIS数据的完整流程
	
	参数:
	hdf_files: HDF文件列表
	output_dir: 输出目录路径
	region_shp: 区域边界Shapefile路径 (可选)
	
	返回:
	final_output: 最终的土地覆盖文件路径str
	"""
	print("=== 开始处理MCD12Q1数据 ===")
	
	# 确保输出目录存在
	os.makedirs(output_dir, exist_ok=True)
	print(f"输出目录: {output_dir}")
	
	# 1. 提取HDF数据
	print("1. 提取HDF数据...")
	extracted_tifs = extract_modis_hdf(hdf_files, output_dir)
	
	# 2. 投影转换和拼接
	print("2. 投影转换和拼接...")
	lc_type1_files = [f for f in extracted_tifs if "LC_Type1" in f]
	mosaic_file = reproject_and_mosaic(lc_type1_files, output_dir)
	
	# 3. 裁剪到研究区域 (如果提供了边界文件)
	if region_shp and os.path.exists(region_shp):
		print("3. 裁剪到研究区域...")
		final_output = clip_to_shp_region(mosaic_file, region_shp, output_dir)
	else:
		print("3. 跳过裁剪步骤 (未提供边界文件)")
		final_output = mosaic_file
	
	print("\n=== 数据处理完成 ===")
	print(f"最终土地覆盖文件: {final_output}")
	
	return final_output

# 完整使用流程示例
if __name__ == "__main__":
	# 参数设置
	hdf_files = glob("MODIS/MCD12Q1*.hdf")  # HDF文件路径模式
	output_dir = "./modis_data"  # 输出目录
	# region_shp = "../grid/Lower_YRB.shp"  # 研究区域边界文件 (可选)
	
	# 执行完整流程
	final_lc_file = process_modis_data(
		hdf_files=hdf_files,
		output_dir=output_dir
	)
	
	# 输出结果摘要
	print(f"\n输出文件:")
	print(f"  - 土地覆盖数据: {final_lc_file}")
	
	# 显示文件信息
	dataset = gdal.Open(final_lc_file)
	print(f"\n文件信息:")
	print(f"  尺寸: {dataset.RasterXSize} x {dataset.RasterYSize}")
	print(f"  波段数: {dataset.RasterCount}")
	print(f"  地理范围: {dataset.GetGeoTransform()}")
	dataset = None