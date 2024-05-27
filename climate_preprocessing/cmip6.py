#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 13:41:04 2022

@author: carmensteinmann
"""
from pathlib import Path
import xarray as xr
import os
import sys
import pathlib

dir_preprocessing = os.path.join(str(pathlib.Path().resolve()))
dir_data = os.path.join(dir_preprocessing, 'data')
sys.path.append(dir_preprocessing)
import historic_metric as fct_his



thresholds = [ ]
lon_bounds_t = [235, 294] #longitudinal values range from 0 to 360
lat_bounds = [24, 50]
geo_area = 'USA'
years = (1980, 2022)
nr_years = 40 
quantile = 0.5
crops = ['maize', 'wheat']

"""Get files"""
import os

# Directory paths
historic_dir = "/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/climate_preprocessing/data/input/CMIP6/historic"
future_dir = "/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/climate_preprocessing/data/input/CMIP6/future"

# Initialize a dictionary to store files by climate model
climate_model_files = {}

# Step 1: Determine climate model used in historic files and store historic file path
for root, _, files in os.walk(historic_dir):
    for filename in files:
        parts = filename.split('_')
        if len(parts) >= 2:
            climate_model = parts[1]
            historic_file = os.path.join(root, filename)
            climate_model_files[climate_model] = {"historic_file": historic_file, "future_files": []}

# Step 2: Find future files for each climate model
for climate_model, data in climate_model_files.items():
    for root, _, files in os.walk(future_dir):
        for filename in files:
            if filename.startswith(f"tasmax_{climate_model}"):
                future_file = os.path.join(root, filename)
                data["future_files"].append(future_file)

#climate_model_files


data_t =  Path(dir_data+'/input/CMIP6/ssp119/tasmax_CanESM5_ssp119_r1i1p1f1_2015-2100.nc') 


t_monthly = fct_his.sm_metrics(data_t, lon_bounds_t, lat_bounds, crops, variable='tasmax')
t_monthly = t_monthly.drop_vars('height')



#mean over the chosen time period
final_year = t_monthly.year.max()
start_year = final_year-nr_years+1
t_sub = t_monthly.sel(year=slice(start_year, final_year))
future_mean = t_sub.mean(dim='year')

# compute the quantile
t_quantile = t_sub.quantile(q=quantile, dim='year') 
binary_data = xr.where(t_sub >= t_quantile, 1, 0)
nr_exceedances = binary_data.sum('year')


"""co-occurence"""
# Initialize an empty dictionary to store the results
concurrence_results = {}

# Loop over each crop to compute concurrences
for crop in crops:
    spring_key = f"{crop}_spring"
    summer_key = f"{crop}_summer"
    
    # Logical AND for spring and summer
    both = (binary_data[spring_key] == 1) & (binary_data[summer_key] == 1)
    
    # Count the number of concurrences per lat/lon by summing along the 'year' dimension
    concurrence = both.sum(dim='year')
    
    # Store the result in the dictionary
    concurrence_results[f"{crop}_concurrence"] = concurrence

# Create a new Dataset to store these results
concurrence_dataset = xr.Dataset(concurrence_results, coords={"lat": binary_data.lat, "lon": binary_data.lon})

# Print the results
print(concurrence_dataset)







#HISTORIC COMPUTATIONS
# import os
# import sys
# import pathlib
# from pathlib import Path
# import xarray as xr
# dir_preprocessing = os.path.join(str(pathlib.Path().resolve()))
# sys.path.append(dir_preprocessing)
# import historic_metric as fct_his

# dir_data = os.path.join(dir_preprocessing, 'data')
# data_t = Path(os.path.join(dir_data, 'input', 'CPC'))
# data_m = Path(os.path.join(dir_data, 'input', 'SMroot_1980-2021_GLEAM_v3.6a_daily_remap05.nc'))
# months = [3, 4, 5, 6, 7, 8]
# thresholds = [15, 20.5, 28.8, 30.5]
# lon_bounds_t = [235, 294] #longitudinal values range from 0 to 360
# lon_bounds_m = [-125, -66] #longitudinal values range from -180 to 180
# lat_bounds = [50, 24]
# geo_area = 'USA'
# years = (1980, 2022)
# monthly_t_metrics(data_t, lon_bounds_t, lat_bounds, months, thresholds, years)