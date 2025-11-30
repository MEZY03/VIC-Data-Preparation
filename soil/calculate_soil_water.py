##########################
# calculate_soil_water.py
##########################
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def calculate_soil_water_characteristic(Sand, Clay, OrganicMatter, Salinity=0, Gravel=0, Compaction=1.0):
	"""
	基于 Saxton & Rawls (2006) 公式计算土壤水力参数
	计算结果应与SPAW软件的输出一致
	
	参数:
	Sand: 砂粒含量，质量百分比 (0.1-99.9)
	Clay: 粘粒含量，质量百分比 (0.1-99.9)
	OrganicMatter: 有机质含量，质量百分比 (0.1-20)
	Salinity: 盐度 (dS/m)，默认值为 0
	Gravel: 砾石含量，体积百分比 (0-80)，默认值为 0
	Compaction: 压实因子 (0.9-1.3)，默认值为 1.0
	
	返回:
	soil_water_dict: 包含计算出的土壤水力参数的字典dict
	"""
	
	# 限制输入范围并转换单位
	S = np.clip(Sand, 0.1, 99.9) / 100.0   # 转换为小数
	C = np.clip(Clay, 0.1, 99.9) / 100.0   # 转换为小数
	OM = np.clip(OrganicMatter, 0.1, 20)   # 保持为百分比
	R_v = np.clip(Gravel, 0, 80) / 100.0   # 转换为小数
	DF = np.clip(Compaction, 0.9, 1.3)     # 压实因子

	# 1. 计算基本水分特征参数 (方程1-5)
	# 萎蔫点 θ_1500 (方程1)
	theta_1500t = (-0.024 * S + 0.487 * C + 0.006 * OM +
				   0.005 * (S * OM) - 0.013 * (C * OM) +
				   0.068 * (S * C) + 0.031)
	theta_1500 = theta_1500t + (0.14 * theta_1500t - 0.02)

	# 田间持水量 θ_33 (方程2)
	theta_33t = (-0.251 * S + 0.195 * C + 0.011 * OM +
				 0.006 * (S * OM) - 0.027 * (C * OM) +
				 0.452 * (S * C) + 0.299)
	theta_33 = theta_33t + (1.283 * theta_33t**2 - 0.374 * theta_33t - 0.015)

	# 饱和-田间持水量差 θ_S-33 (方程3)
	theta_S33t = (0.278 * S + 0.034 * C + 0.022 * OM -
				  0.018 * (S * OM) - 0.027 * (C * OM) -
				  0.584 * (S * C) + 0.078)
	theta_S33 = theta_S33t + (0.636 * theta_S33t - 0.107)

	# 饱和含水量 θ_S (方程5)
	theta_S = theta_33 + theta_S33 - 0.097 * S + 0.043

	# 2. 计算容重和密度效应 (方程6-10)
	# 名义容重 (方程6)
	bulk_density_normal = (1 - theta_S) * 2.65  # g/cm³

	# 校正后的容重 (方程7)
	bulk_density_DF = bulk_density_normal * DF  # g/cm³

	# 校正后的饱和含水量 (方程8)
	theta_S_DF = 1 - (bulk_density_DF / 2.65)  # 小数

	# 校正后的田间持水量 (方程9)
	theta_33_DF = theta_33 - 0.2 * (theta_S - theta_S_DF)  # 小数

	# 3. 计算进气吸力 ψ_e (方程4)
	psi_et = (-21.67 * S - 27.93 * C - 81.97 * theta_S33 +
			  71.12 * (S * theta_S33) + 8.29 * (C * theta_S33) +
			  14.05 * (S * C) + 27.16)
	psi_e = psi_et + (0.02 * psi_et**2 - 0.113 * psi_et - 0.70)  # kPa

	# 4. 计算饱和导水率 K_s (方程16, 15, 18)
	# 计算 B 参数 (方程15)
	B = (np.log(1500) - np.log(33)) / (np.log(theta_33_DF) - np.log(theta_1500))
	lambda_val = 1 / B  # 方程18

	# 饱和导水率 (方程16)
	K_s = 1930 * (theta_S_DF - theta_33_DF)**(3 - lambda_val)  # mm/h

	# 5. 应用砾石校正 (方程19-22)
	alpha = bulk_density_DF / 2.65  # α = ρ_DF / 2.65
	R_w = R_v/(R_v * (1 - alpha) + alpha) # 方程19变形

	# 砾石校正后的容重 (方程20)
	bulk_density_gravel = bulk_density_DF * (1 - R_v) + (R_v * 2.65)  # g/cm³

	# 砾石校正后的导水率 (方程22)
	K_b_Ks = (1 - R_w) / (1 - R_w * (1 - 1.5 * alpha))
	K_b = K_b_Ks * K_s  # mm/h

	# 6. 盐度校正 (方程23-24)
	EC = Salinity  # 电导率 dS/m
	psi_o = 36 * EC  # 方程23 - 饱和时的渗透势 (kPa)

	# 7. 计算有效含水量
	available_water = theta_33_DF - theta_1500  # 小数
	
	soil_water_dict = {
		'wilting_point': np.clip(theta_1500, 0.01, 0.5),           # m³/m³
		'field_capacity': np.clip(theta_33_DF, 0.05, 0.6),         # m³/m³
		'saturation': np.clip(theta_S_DF, 0.3, 0.7),               # m³/m³
		'available_water': np.clip(available_water, 0.01, 0.4),    # m³/m³
		'sat_hydraulic_cond': np.clip(K_b, 0.1, 1000),             # mm/h
		'matric_bulk_density': np.clip(bulk_density_DF, 1.0, 2.0), # g/cm³
		'b_parameter': np.clip(B, 1, 20),                          # 无量纲
		'lambda_param': np.clip(lambda_val, 0.05, 2.0),            # 无量纲
		'air_entry_tension': np.clip(psi_e, 1, 50),                # kPa
		'bulk_density_normal': bulk_density_normal,                # g/cm³
		'compaction_factor': DF                                    # 无量纲 (压实因子)
	}
	return soil_water_dict


def process_soil_water_characteristic(input_file, output_file, default_BD=None, default_Compaction=1.0):
	"""
	处理土壤参数计算的主函数
	
	参数:
	input_file: 输入CSV文件路径
	output_file: 输出CSV文件路径
	default_BD: 默认土壤容重 (g/cm³)，如果HWSD数据中无容重数据，则使用此值
	default_Compaction: 默认压实因子，默认值为 1.0
	
	返回:
	soil_water_df: 包含计算出的土壤水力参数的数据框DataFrame
	"""
	
	print("=== 土壤参数计算 ===")
	
	# 1. 读取CSV文件
	print("1. 读取CSV文件...")
	df = pd.read_csv(input_file)
	print(f"  读取到 {len(df)} 条记录")

	# 2. 检查必要的字段
	print("2. 检查必要字段...")
	required_columns = ['grid_id', 'lon_center', 'lat_center', 'T_SAND', 'T_CLAY', 'T_OC']
	missing_columns = [col for col in required_columns if col not in df.columns]

	if missing_columns:
		print(f"错误：缺少必要的字段: {missing_columns}")
		return None

	# 3. 数据清洗和预处理
	print("3. 数据预处理...")
	df_clean = df.copy()

	# 处理数值型字段
	numeric_columns = ['T_SAND', 'T_CLAY', 'T_OC', 'T_GRAVEL', 'T_ECE', 'T_REF_BULK_DENSITY']
	for col in numeric_columns:
		if col in df_clean.columns:
			df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce')
		else:
			if col == 'T_ECE':  # 盐度，单位 dS/m
				df_clean[col] = 0
			elif col == 'T_GRAVEL':  # 砾石，单位 %w
				df_clean[col] = 0
			elif col == 'T_REF_BULK_DENSITY':  # 容重，单位 g/cm³
				df_clean[col] = default_BD

	# 确保坐标字段为数值型
	df_clean['lon_center'] = pd.to_numeric(df_clean['lon_center'], errors='coerce')
	df_clean['lat_center'] = pd.to_numeric(df_clean['lat_center'], errors='coerce')

	# 计算有机质含量 (%w)
	df_clean['T_OM'] = df_clean['T_OC'] * 1.724

	# 4. 批量计算土壤水力参数
	print("4. 计算土壤水力参数...")
	results = []
	error_count = 0
	
	for idx, row in df_clean.iterrows():
		try:
			# 获取容重
			BD = row.get('T_REF_BULK_DENSITY', default_BD)

			params = calculate_soil_water_characteristic(
				Sand=row['T_SAND'],              # %w
				Clay=row['T_CLAY'],              # %w
				OrganicMatter=row['T_OM'],       # %w
				Salinity=row.get('T_ECE', 0),    # dS/m
				Gravel=row.get('T_GRAVEL', 0),   # %w
				Compaction=default_Compaction,   # 压实因子
			)

			# 组织结果
			result_row = {
				'grid_id': row['grid_id'],
				'lon_center': row['lon_center'],
				'lat_center': row['lat_center'],
				'T_SAND': row['T_SAND'],                         # %w
				'T_CLAY': row['T_CLAY'],                         # %w
				'T_SILT': 100 - row['T_SAND'] - row['T_CLAY'],   # %w
				'T_OC': row['T_OC'],                             # %w
				'T_OM': row['T_OM'],                             # %w
				'T_GRAVEL': row.get('T_GRAVEL', 0),              # %w
				'T_ECE': row.get('T_ECE', 0),                    # dS/m
				'T_REF_BULK_DENSITY': BD                         # g/cm³
			}
			result_row.update(params)
			results.append(result_row)

		except Exception as e:
			error_count += 1
			if error_count <= 5:  # 只打印前5个错误
				print(f"  警告: 计算第 {idx} 行时出错: {e}")
			continue

	if error_count > 0:
		print(f"  共有 {error_count} 个网格计算失败")

	# 5. 创建结果DataFrame并保存结果
	print("5. 保存结果...")
	soil_water_df = pd.DataFrame(results)
	soil_water_df.to_csv(output_file, index=False)
	print(f"  土壤参数文件: {output_file}")
	print(f"  成功计算了 {len(soil_water_df)} 个网格的土壤参数")

	# 6. 输出统计信息
	print("\n=== 土壤参数统计 ===")
	key_parameters = ['wilting_point', 'field_capacity', 'saturation', 'available_water', 'sat_hydraulic_cond']
	for col in key_parameters:
		if col in soil_water_df.columns:
			valid_data = soil_water_df[col].dropna()
			if len(valid_data) > 0:
				if col == 'sat_hydraulic_cond':
					unit = "mm/h"
				else:
					unit = "m³/m³"
				print(f"  {col}: {valid_data.mean():.4f} {unit} ± {valid_data.std():.4f}")

	return soil_water_df


# 完整使用流程示例
if __name__ == "__main__":
	# 参数设置
	input_file = "./soil_data/grid_soil.csv"  # 输入土壤数据文件
	output_file = "./soil_data/grid_soil_water.csv"  # 输出土壤水力文件
	
	# 执行土壤参数计算
	soil_water_df = process_soil_water_characteristic(input_file, output_file)
	
	# 输出结果摘要
	print(f"\n输出文件:")
	print(f"  - 土壤参数: {output_file}")
	
	# 显示前几个网格的土壤水力参数
	print(f"\n前5个网格的土壤水力参数:")
	display_columns = ['grid_id', 'lon_center', 'lat_center', 
					  'wilting_point', 'field_capacity', 'saturation', 'sat_hydraulic_cond']
	available_columns = [col for col in display_columns if col in soil_water_df.columns]
	print(soil_water_df[available_columns].head(5).to_string(index=False))