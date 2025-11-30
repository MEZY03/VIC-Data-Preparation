##########################
# create_forcing.py
##########################
import xarray as xr
import numpy as np
import pandas as pd
from scipy import interpolate
import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

class VICForcingGenerator:
    """
    VIC模型驱动数据生成器
    """
    
    def __init__(self, config: Dict):
        """
        初始化驱动数据生成器
        
        参数:
        config: 配置字典, 包含文件路径和参数设置
        """
        self.config = config
        self.vic_variables = {}
        self.grid_info = {}
        self.datasets = {}
        
    @staticmethod
    def calculate_vapor_pressure(specific_humidity: np.ndarray, surface_pressure: np.ndarray) -> np.ndarray:
        """
        根据比湿和地表气压计算水汽压
        """
        vapor_pressure_pa = specific_humidity * surface_pressure / (0.622 + 0.378 * specific_humidity)
        vapor_pressure_kpa = vapor_pressure_pa / 1000.0
        return vapor_pressure_kpa
    
    @staticmethod
    def estimate_incoming_longwave_radiation(air_temperature: np.ndarray) -> np.ndarray:
        """
        简化版长波辐射计算 - 使用固定发射率0.85
        
        参数:
        air_temperature: 气温 (K)
        
        返回:
        longwave_radiation: 下行长波辐射通量 (W/m²)
        """
        stefan_boltzmann = 5.67e-8  # Stefan-Boltzmann常数
        emissivity = 0.85  # 固定发射率
        
        longwave_radiation = emissivity * stefan_boltzmann * air_temperature**4
        return longwave_radiation
    
    def interpolate_to_target_grid(self, source_data: xr.DataArray, method: str = 'linear') -> np.ndarray:
        """
        将气象数据插值到目标网格系统
        """
        source_lons = source_data.lon.values
        source_lats = source_data.lat.values
        
        n_time_steps = len(source_data.time)
        n_grid_cells = len(self.grid_info['longitudes'])
        interpolated_data = np.full((n_time_steps, n_grid_cells), np.nan)
        
        for time_index in range(n_time_steps):
            interpolation_function = interpolate.RegularGridInterpolator(
                (source_lats, source_lons),
                source_data[time_index, :, :].values,
                method=method,
                bounds_error=False,
                fill_value=np.nan
            )
            
            target_coordinates = np.column_stack((self.grid_info['latitudes'], self.grid_info['longitudes']))
            interpolated_data[time_index, :] = interpolation_function(target_coordinates)
        return interpolated_data
    
    def load_grid_information(self) -> bool:
        """
        加载目标网格信息
        """
        print("1. 加载目标网格信息...")
        
        try:
            soil_params = pd.read_csv(self.config['soil_param_file'], delim_whitespace=True, header=None)
            
            self.grid_info['ids'] = soil_params[1].values
            self.grid_info['latitudes'] = soil_params[2].values
            self.grid_info['longitudes'] = soil_params[3].values
            
            print(f"  网格数量: {len(self.grid_info['latitudes'])}")
            print(f"  纬度范围: {np.min(self.grid_info['latitudes']):.3f} ~ {np.max(self.grid_info['latitudes']):.3f}")
            print(f"  经度范围: {np.min(self.grid_info['longitudes']):.3f} ~ {np.max(self.grid_info['longitudes']):.3f}")
            return True
            
        except Exception as e:
            print(f"  错误: 无法加载网格信息 - {e}")
            return False
    
    def load_meteorological_data(self) -> bool:
        """
        加载气象强迫数据
        """
        print("2. 加载气象强迫数据...")
        
        try:
            # 加载必需的气象变量
            required_files = [
                ('temperature', 'temperature_file'),
                ('pressure', 'pressure_file'),
                ('humidity', 'humidity_file'),
                ('wind', 'wind_file'),
                ('shortwave', 'shortwave_file'),
                ('precipitation', 'precipitation_file')
            ]
            
            for ds_name, config_key in required_files:
                if config_key in self.config and self.config[config_key]:
                    self.datasets[ds_name] = xr.open_dataset(self.config[config_key])
                    print(f"  {ds_name}: {os.path.basename(self.config[config_key])}")
            
            # 可选的长波辐射
            if 'longwave_file' in self.config and self.config['longwave_file']:
                self.datasets['longwave'] = xr.open_dataset(self.config['longwave_file'])
                print(f"  longwave: {os.path.basename(self.config['longwave_file'])}")
            else:
                print("  longwave: 使用简化公式估算")
            
            # 获取时间信息
            if self.datasets:
                first_ds = next(iter(self.datasets.values()))
                self.time_coordinates = first_ds.time.values
                self.n_time_steps = len(self.time_coordinates)
                
                print(f"  时间步长: {self.n_time_steps}")
                print(f"  时间范围: {self.time_coordinates[0]} ~ {self.time_coordinates[-1]}")
                return True
            else:
                print("  错误: 未加载任何气象数据")
                return False
                
        except Exception as e:
            print(f"  错误: 加载气象数据失败 - {e}")
            return False
    
    def _apply_unit_conversion(self, data: np.ndarray, conversion_key: str) -> np.ndarray:
        """
        应用单位转换函数
        
        参数:
        data: 输入数据
        conversion_key: 单位转换配置键名
        
        返回:
        转换后的数据
        """
        if conversion_key in self.config and self.config[conversion_key] is not None:
            try:
                return self.config[conversion_key](data)
            except Exception as e:
                print(f"  警告: 单位转换失败 {conversion_key} - {e}")
                return data
        return data
    
    def process_temperature_data(self) -> bool:
        """处理气温数据"""
        if 'temperature' not in self.datasets:
            return False
            
        print("  处理气温数据...")
        temperature_data = self.datasets['temperature'][self.config['temperature_var']]
        temperature_interp = self.interpolate_to_target_grid(temperature_data)
        
        # 应用单位转换
        temperature_interp = self._apply_unit_conversion(temperature_interp, 'temperature_conv')
        self.vic_variables['AIR_TEMP'] = temperature_interp
        return True
    
    def process_precipitation_data(self) -> bool:
        """处理降水数据"""
        if 'precipitation' not in self.datasets:
            return False
            
        print("  处理降水数据...")
        precipitation_data = self.datasets['precipitation'][self.config['precipitation_var']]
        precipitation_interp = self.interpolate_to_target_grid(precipitation_data)
        
        # 应用单位转换
        precipitation_interp = self._apply_unit_conversion(precipitation_interp, 'precipitation_conv')
        self.vic_variables['PREC'] = precipitation_interp
        return True
    
    def process_pressure_data(self) -> bool:
        """处理气压数据"""
        if 'pressure' not in self.datasets:
            return False
            
        print("  处理气压数据...")
        pressure_data = self.datasets['pressure'][self.config['pressure_var']]
        pressure_interp = self.interpolate_to_target_grid(pressure_data)
        
        # 应用单位转换
        pressure_interp = self._apply_unit_conversion(pressure_interp, 'pressure_conv')
        self.vic_variables['PRESSURE'] = pressure_interp
        return True
    
    def process_shortwave_radiation(self) -> bool:
        """处理短波辐射数据"""
        if 'shortwave' not in self.datasets:
            return False
            
        print("  处理短波辐射...")
        shortwave_data = self.datasets['shortwave'][self.config['shortwave_var']]
        shortwave_interp = self.interpolate_to_target_grid(shortwave_data)
        
        # 应用单位转换
        shortwave_interp = self._apply_unit_conversion(shortwave_interp, 'shortwave_conv')
        self.vic_variables['SWDOWN'] = shortwave_interp
        return True
    
    def process_longwave_radiation(self) -> bool:
        """处理长波辐射数据"""
        # 优先使用观测的长波辐射
        if ('longwave' in self.datasets and 
            self.config.get('longwave_var') and 
            self.config['longwave_var'] in self.datasets['longwave'].variables):
            
            print("  处理长波辐射(观测)...")
            longwave_data = self.datasets['longwave'][self.config['longwave_var']]
            longwave_interp = self.interpolate_to_target_grid(longwave_data)
            
            # 应用单位转换
            longwave_interp = self._apply_unit_conversion(longwave_interp, 'longwave_conv')
            self.vic_variables['LWDOWN'] = longwave_interp
            return True
        else:
            # 使用简化公式估算长波辐射
            print("  使用简化公式估算长波辐射...")
            if 'temperature' in self.datasets:
                temperature_data = self.datasets['temperature'][self.config['temperature_var']]
                temperature_k_interp = self.interpolate_to_target_grid(temperature_data)
                
                # 使用简化公式计算长波辐射
                estimated_longwave = self.estimate_incoming_longwave_radiation(temperature_k_interp)
                self.vic_variables['LWDOWN'] = estimated_longwave
                return True
        
        print("  警告: 无法处理长波辐射数据")
        return False
    
    def process_humidity_data(self) -> bool:
        """处理湿度数据并计算水汽压"""
        if 'humidity' not in self.datasets or 'PRESSURE' not in self.vic_variables:
            return False
            
        print("  计算水汽压...")
        humidity_data = self.datasets['humidity'][self.config['humidity_var']]
        humidity_interp = self.interpolate_to_target_grid(humidity_data)
        
        # 应用单位转换
        humidity_interp = self._apply_unit_conversion(humidity_interp, 'humidity_conv')
        
        # 计算水汽压 (注意单位: PRESSURE是kPa, 需要转换为Pa)
        vapor_pressure = self.calculate_vapor_pressure(
            humidity_interp, 
            self.vic_variables['PRESSURE'] * 1000
        )
        self.vic_variables['VP'] = vapor_pressure
        return True
    
    def process_wind_data(self) -> bool:
        """处理风速数据"""
        if 'wind' not in self.datasets:
            return False
            
        print("  处理风速数据...")
        wind_data = self.datasets['wind'][self.config['wind_var']]
        wind_interp = self.interpolate_to_target_grid(wind_data)
        
        # 应用单位转换
        wind_interp = self._apply_unit_conversion(wind_interp, 'wind_conv')
        self.vic_variables['WIND'] = wind_interp
        return True
    
    def generate_vic_forcing_files(self) -> bool:
        """
        生成VIC格式的驱动文件
        """
        print("4. 生成VIC驱动文件...")
        
        output_dir = self.config['output_dir']
        os.makedirs(output_dir, exist_ok=True)
        
        # VIC变量顺序 (固定)
        vic_var_order = ['AIR_TEMP', 'PREC', 'PRESSURE', 'SWDOWN', 'LWDOWN', 'VP', 'WIND']
        
        n_grid_cells = len(self.grid_info['latitudes'])
        n_time_steps = self.n_time_steps
        
        print(f"  输出目录: {output_dir}")
        print(f"  网格点数: {n_grid_cells}")
        print(f"  时间步数: {n_time_steps}")
        
        # 为每个网格点生成驱动文件
        file_count = 0
        for grid_index in range(n_grid_cells):
            current_lat = self.grid_info['latitudes'][grid_index]
            current_lon = self.grid_info['longitudes'][grid_index]
            current_grid_id = self.grid_info['ids'][grid_index]
            
            # 构建文件名
            filename = f"forcing_{current_lat:.6f}_{current_lon:.6f}"
            file_path = os.path.join(output_dir, filename)
            
            try:
                # 使用二进制模式写入, 避免编码问题
                with open(file_path, 'w', encoding='utf-8') as file_handler:
                    for time_index in range(n_time_steps):
                        data_line = []
                        for variable_name in vic_var_order:
                            if variable_name in self.vic_variables:
                                variable_data = self.vic_variables[variable_name]
                                value = float(variable_data[time_index, grid_index])
                                
                                # 数据有效性检查
                                if np.isnan(value) or np.isinf(value):
                                    value = 0.0
                                
                                data_line.append(f"{value:.4f}")
                            else:
                                data_line.append("0.0000")
                        
                        file_handler.write(" ".join(data_line) + "\n")
                
                file_count += 1
                if file_count % 100 == 0:
                    print(f"  已生成 {file_count}/{n_grid_cells} 个文件")
                    
            except Exception as e:
                print(f"  警告: 生成文件 {filename} 失败 - {e}")
                continue
        
        # 生成网格信息元数据文件
        self._generate_grid_metadata_file(output_dir)
        
        print(f"  共生成 {file_count}/{n_grid_cells} 个驱动文件")
        return file_count > 0
    
    def _generate_grid_metadata_file(self, output_dir: str):
        """生成网格信息元数据文件"""
        metadata_filepath = os.path.join(output_dir, "vic_forcing_info.txt")
        
        try:
            with open(metadata_filepath, 'w', encoding='utf-8') as meta_file:
                # 写入文件头信息
                meta_file.write("# VIC模型驱动数据网格信息文件\n")
                meta_file.write(f"# 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                meta_file.write("# 列说明: 网格ID 纬度(度) 经度(度) 文件名\n")
                meta_file.write("# 变量顺序: AIR_TEMP(°C) PREC(mm) PRESSURE(kPa) SWDOWN(W/m²) LWDOWN(W/m²) VP(kPa) WIND(m/s)\n\n")
                
                # 写入每个网格的信息
                for i in range(len(self.grid_info['latitudes'])):
                    filename = f"forcing_{self.grid_info['latitudes'][i]:.6f}_{self.grid_info['longitudes'][i]:.6f}"
                    meta_file.write(f"{self.grid_info['ids'][i]} {self.grid_info['latitudes'][i]:.8f} {self.grid_info['longitudes'][i]:.8f} {filename}\n")
            
            print(f"  网格元数据文件: {metadata_filepath}")
        except Exception as e:
            print(f"  警告: 生成网格元数据文件失败 - {e}")
    
    def close_datasets(self):
        """关闭所有打开的数据集"""
        for name, dataset in self.datasets.items():
            try:
                dataset.close()
            except:
                pass
    
    def generate_forcing_data(self) -> bool:
        """
        生成VIC驱动数据的主函数
        """
        print("=== VIC模型驱动数据生成 ===")
        start_time = datetime.now()
        
        try:
            # 1. 加载网格信息
            if not self.load_grid_information():
                return False
            
            # 2. 加载气象数据
            if not self.load_meteorological_data():
                return False
            
            # 3. 处理各个气象变量
            print("3. 处理气象变量...")
            
            # 处理顺序很重要：气压和温度需要先处理, 因为其他变量依赖它们
            processing_steps = [
                (self.process_temperature_data, "气温"),
                (self.process_pressure_data, "气压"),
                (self.process_precipitation_data, "降水"),
                (self.process_shortwave_radiation, "短波辐射"),
                (self.process_humidity_data, "湿度"),
                (self.process_wind_data, "风速"),
                (self.process_longwave_radiation, "长波辐射")
            ]
            
            success_count = 0
            for process_func, desc in processing_steps:
                if process_func():
                    success_count += 1
                    print(f"  {desc}处理完成")
                else:
                    print(f"  {desc}处理跳过")
            
            # 检查必需变量
            essential_vars = ['AIR_TEMP', 'PREC', 'PRESSURE', 'SWDOWN']
            missing_essential = [var for var in essential_vars if var not in self.vic_variables]
            
            if missing_essential:
                print(f"  错误: 缺少必需变量: {missing_essential}")
                return False
            
            # 4. 生成驱动文件
            if not self.generate_vic_forcing_files():
                return False
            
            # 5. 清理资源
            self.close_datasets()
            
            processing_time = datetime.now() - start_time
            print(f"=== 驱动数据生成完成! 总耗时: {processing_time} ===")
            return True
            
        except Exception as e:
            print(f"=== 驱动数据生成失败: {e} ===")
            import traceback
            traceback.print_exc()
            self.close_datasets()
            return False

# 完整使用流程示例
if __name__ == "__main__":
    print("VIC模型驱动数据生成工具")
    print("功能: 将气象强迫数据转换为VIC模型输入格式")
    
    # 创建配置
    config = {
        # 目标网格信息
        'soil_param_file': "../soil/soil_param.txt",
        
        # 输入文件配置
        'temperature_file': "./forcing_data/temperature.nc",
        'precipitation_file': "./forcing_data/precipitation.nc",
        'pressure_file': "./forcing_data/pressure.nc",
        'shortwave_file': "./forcing_data/shortwave_radiation.nc",
        'longwave_file': "./forcing_data/longwave_radiation.nc",  # 可选
        'humidity_file': "./forcing_data/specific_humidity.nc", 
        'wind_file': "./forcing_data/wind_speed.nc",
        
        # 变量名映射
        'temperature_var': 'TLML',
        'precipitation_var': 'PRECTOT',
        'pressure_var': 'PS',
        'shortwave_var': 'SWGDN',
        'longwave_var': 'LWGAB',  # 可选
        'humidity_var': 'QLML',
        'wind_var': 'SPEED',
        
        # 单位转换函数 (如果前处理已完成转换, 设为None)
        'temperature_conv': lambda x:x-273.15,    # K -> °C
        'precipitation_conv': lambda x:x*3600,    # kg/m²/s -> mm/h
        'pressure_conv': lambda x:x/1000,         # Pa -> kPa  
        'shortwave_conv': None,        # 通常无需转换
        'longwave_conv': None,         # 通常无需转换
        'humidity_conv': None,         # 通常无需转换
        'wind_conv': None,             # 通常无需转换
        
        # 输出配置
        'output_dir': "./forcings"
    }
    
    # 初始化生成器
    generator = VICForcingGenerator(config)
    
    # 生成驱动数据
    success = generator.generate_forcing_data()
    
    # 输出结果
    if success:
        print(f"\n生成成功, 输出文件:")
        print(f"  - 驱动数据文件: {config['output_dir']}/forcing_*")
        print(f"  - 网格信息文件: {config['output_dir']}/vic_forcing_info.txt")
        print(f"\n下一步: 使用生成的强迫数据和配置文件运行VIC模型")
    else:
        print("\n生成失败, 请检查配置和输入数据")