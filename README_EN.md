## Project Introduction
This repository provides a set of Python script tools for preparing input data for the VIC (Variable Infiltration Capacity) hydrological model. The project aims to automate the entire process from raw data downloading to generating VIC-recognized parameter files, helping researchers quickly set up VIC model experimental environments.

For complete detailed tutorials, data processing principles, and important notes about the project, please refer to my CSDN blog post: [Complete Tutorial on VIC Model Input Data Preparation](https://).

## Key Preparation Steps
1. Grid Creation: Generate simulation grid systems based on the study area Shapefile (grid_coordinates.csv).

2. Topographic Data: Extract grid elevation from DEM data (e.g., Copernicus GLO-30).

3. Climatic Precipitation: Process long-term time series precipitation data (e.g., MERRA-2) to generate climatological annual precipitation background fields.

4. Soil Parameters: Calculate soil hydraulic parameters using Saxton & Rawls (2006) formulas based on HWSD soil data, and integrate topographic and precipitation information to generate soil_param.txt.

5. Vegetation Parameters: Calculate vegetation type proportions based on MODIS land cover data (MCD12Q1) and generate veg_param.txt.

6. Meteorological Forcing Data: Interpolate hourly meteorological data (e.g., MERRA-2) to grid points, generating forcing_[lat]_[lon].txt files for each grid.

7. One-click Run Script:
```
./run_all_scripts.sh
```

## Dependencies
Core Python Libraries: Python 3.11+ Standard Libraries

You can install the main Python dependencies using the following command:
```
conda install numpy pandas scipy xarray matplotlib gdal libgdal-hdf4 geopandas shapely rasterstats netcdf4 h5netcdf
```
Auxiliary Tools: CDO, NCO

## Official VIC Resources
[VIC Official Documentation](https://vic.readthedocs.io/en/master/)

[VIC Official GitHub Repository](https://github.com/UW-Hydro/VIC)

[VIC Sample Data](https://github.com/UW-Hydro/VIC_sample_data)
