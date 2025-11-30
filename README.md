## 项目简介

本仓库提供了一套用于制备 **VIC (Variable Infiltration Capacity)** 水文模型输入数据的Python脚本工具集。项目旨在自动化处理从原始数据下载到生成VIC可识别参数文件的全流程，帮助研究人员快速构建VIC模型实验环境。

关于项目的**完整详细教程**、数据处理原理及注意事项，请参阅我在CSDN上发布的博文：[《VIC模型输入数据制备完整教程》](https://blog.csdn.net/your_blog_article_link) 。

## 关键制备步骤

1. 网格创建：基于研究区域Shapefile生成模拟网格系统 (grid_coordinates.csv)。

2. 地形数据：从DEM数据（如Copernicus GLO-30）中提取网格高程。

3. 气候态降水：处理长时间序列降水数据（如MERRA-2）以生成气候态背景场。

4. 土壤参数：基于HWSD土壤数据，应用Saxton & Rawls (2006)公式计算土壤水力学参数，并整合地形与降水信息生成 soil_param.txt。

5. 植被参数：基于MODIS土地覆盖数据 (MCD12Q1) 计算植被类型比例，生成 veg_param.txt。

6. 气象强迫数据：将逐小时气象数据（如MERRA-2）插值到网格点，生成各网格的 forcing_[lat]_[lon].txt 文件。

7. 一键运行脚本：
```
./run_all_scripts.sh
```

## 依赖环境
核心Python库：Python 3.11+标准库

您可以通过以下命令安装主要Python依赖项：
```
conda install numpy pandas scipy xarray matplotlib gdal libgdal-hdf4 geopandas shapely rasterstats netcdf4 h5netcdf
```

辅组工具：CDO、NCO

## 官方VIC资源
[VIC官方文档](https://github.com/UW-Hydro/VIC)

[VIC官方GitHub仓库](https://github.com/UW-Hydro/VIC)

[VIC示例数据](https://github.com/UW-Hydro/VIC_sample_data)
