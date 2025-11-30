##########################
# create_soil_param.py
##########################
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def generate_soil_param(soil_file, dem_file, precip_file, output_soil_file, soil_depth = [0.1, 0.3, 1.5]):
	"""
	生成VIC模型所需的土壤参数文件
	
	参数:
	soil_file: 土壤参数CSV文件路径
	dem_file: 地形数据CSV文件路径  
	precip_file: 降水数据CSV文件路径
	output_soil_file: 输出的土壤参数文件路径
	soil_depth: 各层土壤厚度，默认为3层
	"""
	print("=== 生成VIC土壤参数文件 ===")
	
	# 1. 读取数据
	print("1. 读取数据文件...")
	soil_df = pd.read_csv(soil_file)
	dem_df = pd.read_csv(dem_file)
	precip_df = pd.read_csv(precip_file)
	
	print(f"  土壤数据: {len(soil_df)} 个网格")
	print(f"  地形数据: {len(dem_df)} 个网格") 
	print(f"  降水数据: {len(precip_df)} 个网格")
	
	# 2. 合并数据
	print("2. 合并数据...")
	merged_df = pd.merge(soil_df, dem_df, on='grid_id', how='inner')
	merged_df = pd.merge(merged_df, precip_df, on='grid_id', how='inner')
	print(f"  合并后: {len(merged_df)} 个网格")
	
	# 3. 土壤数据清洗
	print("3. 土壤数据清洗...")
	initial_count = len(merged_df)
	
	# 移除缺失值
	merged_df = merged_df.dropna(subset=['T_SAND', 'T_CLAY', 'T_OC'])
	
	# 移除异常值
	merged_df = merged_df[
		(merged_df['T_SAND'] >= 0) & (merged_df['T_SAND'] <= 100) &
		(merged_df['T_CLAY'] >= 0) & (merged_df['T_CLAY'] <= 100) & 
		(merged_df['T_OC'] >= 0) & (merged_df['T_OC'] <= 100)
	]
	
	print(f"  清洗后: {len(merged_df)} 个网格 (过滤了 {initial_count - len(merged_df)} 个)")
	
	if len(merged_df) == 0:
		print("错误: 清洗后无有效数据")
		return None
	
	# 4. 填充缺失值
	print("4. 填充缺失值...")
	soil_vic_df = merged_df.copy()
	
	# 为可能缺失的列设置默认值
	default_values = {
		'sat_hydraulic_cond': 10,
		'lambda_param': 0.5, 
		'field_capacity': 0.3,
		'wilting_point': 0.1,
		'saturation': 0.5,
		'air_entry_tension': 10,
		'matric_bulk_density': 1.5,
		'elevation': 100,
		'precipitation': 1000
	}
	
	for col, default in default_values.items():
		if col in soil_vic_df.columns:
			soil_vic_df[col] = soil_vic_df[col].fillna(default)
	
	# 5. 生成土壤参数文件
	print("5. 生成土壤参数文件...")
	
	with open(output_soil_file, 'w') as f:
		for idx, row in soil_vic_df.iterrows():
			# 基本网格信息
			run_cell = 1
			gridcel = int(row['grid_id'])
			lat = row['lat_center']
			lon = row['lon_center']
			
			# 基本水力参数
			wilting_point = row['wilting_point']
			field_capacity = row['field_capacity']
			saturation = row['saturation']
			sat_hydraulic_cond = row['sat_hydraulic_cond']
			lambda_param = row['lambda_param']
			
			# 下渗与基流参数
			infilt = 0.3
			Ds = 0.001
			Ws = 0.8
			c = 2.0
			
			# 土壤水力参数
			expt = 3 + 2 / lambda_param
			Ksat = sat_hydraulic_cond * 24
			slope = 0.1
			Dsmax = Ksat * slope
			phi_s = -99
			
			# 土壤深度
			depth = np.array(soil_depth)
			
			# 初始土壤水分
			init_moist = field_capacity * depth * 1000 * 0.6
			
			# 地形和气候参数
			elev = row['elevation']
			avg_T = 15 - (elev / 100) * 0.6
			dp = 4.0
			
			# 土壤物理参数
			bubble = row['air_entry_tension'] * 10.197
			# bubble = 0.32 * expt + 4.3       # 备选方法
			quartz = 0.4 + 0.003 * row['T_SAND'] 
			bulk_density = row['matric_bulk_density'] * 1000
			soil_density = 2685
			
			# 时区和水分特征参数
			off_gmt = 8
			Wcr_FRACT = 0.7 * field_capacity / saturation
			Wpwp_FRACT = min(wilting_point / saturation, Wcr_FRACT)
			
			# 地表参数
			rough = 0.01
			snow_rough = 0.001
			annual_prec = row['precipitation']
			
			# 其他参数
			resid_moist = wilting_point * 0.3
			fs_active = 0
			
			# 组织输出行 - 按参数类别分组，避免过长行
			# 1. 基本网格信息 (4个参数)
			grid_info = f"{run_cell} {gridcel} {lat:.6f} {lon:.6f}"
			
			# 2. 下渗与基流参数 (5个参数)
			infiltration_params = f"{infilt:.4f} {Ds:.4f} {Dsmax:.4f} {Ws:.4f} {c:.0f}"
			
			# 3. 土壤水力参数 (12个参数)
			hydraulic_params = f"{expt:.4f} {expt:.4f} {expt:.4f} {Ksat:.4f} {Ksat:.4f} {Ksat:.4f}"
			hydraulic_params += f" {phi_s:.0f} {phi_s:.0f} {phi_s:.0f}"
			hydraulic_params += f" {init_moist[0]:.4f} {init_moist[1]:.4f} {init_moist[2]:.4f}"
			
			# 4. 地形和气候参数 (6个参数)
			terrain_params = f"{elev:.2f} {depth[0]:.4f} {depth[1]:.4f} {depth[2]:.4f} {avg_T:.2f} {dp:.0f}"
			
			# 5. 土壤物理参数 (12个参数)
			soil_physics = f"{bubble:.4f} {bubble:.4f} {bubble:.4f}"
			soil_physics += f" {quartz:.3f} {quartz:.3f} {quartz:.3f}"
			soil_physics += f" {bulk_density:.1f} {bulk_density:.1f} {bulk_density:.1f}"
			soil_physics += f" {soil_density:.0f} {soil_density:.0f} {soil_density:.0f}"
			
			# 6. 时区和水分特征参数 (7个参数)
			moisture_params = f"{off_gmt:.0f}"
			moisture_params += f" {Wcr_FRACT:.4f} {Wcr_FRACT:.4f} {Wcr_FRACT:.4f}"
			moisture_params += f" {Wpwp_FRACT:.4f} {Wpwp_FRACT:.4f} {Wpwp_FRACT:.4f}"
			
			# 7. 地表和其他参数 (7个参数)
			surface_params = f"{rough:.4f} {snow_rough:.4f} {annual_prec:.2f}"
			surface_params += f" {resid_moist:.4f} {resid_moist:.4f} {resid_moist:.4f}"
			surface_params += f" {fs_active:.0f}"
			
			# 合并所有参数
			output_line = f"{grid_info} {infiltration_params} {hydraulic_params} {terrain_params} {soil_physics} {moisture_params} {surface_params}"
			
			f.write(output_line + '\n')
	
	print(f"VIC土壤参数文件生成完成: {output_soil_file}")
	print(f"共处理了 {len(soil_vic_df)} 个网格")
	
	# 显示统计信息
	print("\n=== 统计信息 ===")
	print(f"网格数量: {len(soil_vic_df)}")
	print(f"纬度范围: {soil_vic_df['lat_center'].min():.4f} - {soil_vic_df['lat_center'].max():.4f}")
	print(f"经度范围: {soil_vic_df['lon_center'].min():.4f} - {soil_vic_df['lon_center'].max():.4f}")
	print(f"高程范围: {soil_vic_df['elevation'].min():.1f} - {soil_vic_df['elevation'].max():.1f} m")
	print(f"年降水范围: {soil_vic_df['precipitation'].min():.1f} - {soil_vic_df['precipitation'].max():.1f} mm")

def create_soil_param_description(output_desc_file="soil_param_description.txt", soil_depth = [0.1, 0.3, 1.5]):
	"""
	创建土壤参数文件的描述文档
	"""
	description = """# VIC土壤参数文件说明
生成时间: {time}

## 文件格式
总列数: 53列

### 1. 基本网格信息 (4列)
列号\t变量名\t单位\t描述
1\trun_cell\t-\t1=运行网格, 0=不运行
2\tgridcel\t-\t网格编号
3\tlat\t度\t纬度
4\tlon\t度\t经度

### 2. 下渗与基流参数 (5列)
5\tinfilt\t-\t变量下渗曲线参数
6\tDs\t分数\t非线性基流开始的比例
7\tDsmax\tmm/天\t最大基流速度
8\tWs\t分数\t非线性基流发生的土壤水分比例
9\tc\t-\t基流曲线指数

### 3. 土壤水力参数 (12列)
10-12\texpt\t-\tCampbell水力传导指数 (3层)
13-15\tKsat\tmm/天\t饱和水力传导率 (3层)
16-18\tphi_s\t-\t土壤水分扩散参数(未使用) (3层)
19-21\tinit_moist\tmm\t初始土壤水分 (3层)

### 4. 地形和气候参数 (6列)
22\telev\tm\t平均高程
23-25\tdepth\tm\t土壤层厚度 (3层)
26\tavg_T\t°C\t平均土壤温度
27\tdp\tm\t土壤热阻尼深度

### 5. 土壤物理参数 (12列)
28-30\tbubble\tcm\t气泡压力 (3层)
31-33\tquartz\t-\t石英含量 (3层)
34-36\tbulk_density\tkg/m³\t土壤容重 (3层)
37-39\tsoil_density\tkg/m³\t土壤颗粒密度 (3层)

### 6. 时区和水分特征参数 (7列)
40\toff_gmt\t小时\t时区偏移
41-43\tWcr_FRACT\t分数\t临界点水分比例 (3层)
44-46\tWpwp_FRACT\t分数\t萎蔫点水分比例 (3层)

### 7. 地表和其他参数 (7列)
47\trough\tm\t地表粗糙度
48\tsnow_rough\tm\t雪面粗糙度
49\tannual_prec\tmm\t年平均降水
50-52\tresid_moist\t分数\t残余水分含量 (3层)
53\tfs_active\t-\t冻土算法激活标志

## 参数估算说明
- infilt, Ds, Ws, c: 使用典型文献值
- Dsmax: 基于饱和导水率估算 (Ksat * slope)
- avg_T: 基于高程估算(每100m降温0.6°C)
- quartz: 使用典型值0.6
- soil_density: 使用典型值2685 kg/m³
- Wcr_FRACT: 使用典型值0.7
- Wpwp_FRACT: 萎蔫点/饱和含水量
- 粗糙度参数: 使用典型文献值
- 初始土壤水分: 设为田间持水量的80%
- 残余水分含量: 设为萎蔫点的50%

## 土壤层设置
- 第1层: 0-{layer1:.1f}m (厚度: {thickness1:.1f}m)
- 第2层: {layer1:.1f}-{layer2:.1f}m (厚度: {thickness2:.1f}m)  
- 第3层: {layer2:.1f}-{layer3:.1f}m (厚度: {thickness3:.1f}m)
总厚度: {layer3:.1f}m
""".format(time=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
	layer1=soil_depth[0],
	layer2=soil_depth[0] + soil_depth[1],
	layer3=soil_depth[0] + soil_depth[1] + soil_depth[2],
	thickness1=soil_depth[0],
	thickness2=soil_depth[1],
	thickness3=soil_depth[2]
)
	
	with open(output_desc_file, 'w', encoding='utf-8') as f:
		f.write(description)
	
	print(f"参数描述文件已生成: {output_desc_file}")

# 完整使用流程示例
if __name__ == "__main__":
	# 参数设置
	soil_file = "./soil_data/grid_soil_water.csv"  # 输出土壤水力文件
	dem_file = "../elevation/elevation_data/grid_elevation.csv"    # 输入地形数据文件
	precip_file = "../precipitation/precipitation_data/grid_precipitation.csv"  # 输入降水数据文件
	output_soil_file = "soil_param.txt"     # 输出土壤参数文件
	soil_depth = [0.1, 0.3, 1.5]    # 各层土壤厚度
	
	# 生成VIC土壤参数文件
	generate_soil_param(
		soil_file=soil_file,
		dem_file=dem_file,
		precip_file=precip_file,
		output_soil_file=output_soil_file,
		soil_depth = soil_depth
	)
	
	# 输出结果摘要
	print(f"\n输出文件:")
	print(f"  - 土壤参数文件: {output_soil_file}")
	
	# 生成输出参数描述文件
	create_soil_param_description(soil_depth = soil_depth)
	
	# 显示文件前几行作为示例
	print(f"\n前3个网格的土壤参数:")
	with open(output_soil_file, 'r') as f:
		for i, line in enumerate(f):
			if i < 3:
				print(f"第{i+1}行: {line.strip()}")
			else:
				break