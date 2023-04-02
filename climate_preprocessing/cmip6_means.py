#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 10:45:18 2023

@author: carmenst
"""

import os
import xarray as xr
import numpy as np
from glob import iglob
from pathlib import Path
import pandas as pd
import geopandas as gpd
from rasterio import features
from affine import Affine
import rasterstats as rstats
import rasterio as rio


def create_mean(rootdir_glob, nr_years, old_time_period):
    # This will return absolute paths
    file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]
    #output_path = rootdir_glob.replace("input/**/*", "output")
    
    

    for file in file_list:
        new_filename = file.replace("input", "output")
        new_filename = new_filename.replace(old_time_period, "mean")
        
        ds_tmax= xr.open_dataset(file)
        tmax = ds_tmax.tasmax.values[-nr_years:, :, :]
        mean_tmax = np.mean(tmax, axis=0)
        

        ds_output = xr.Dataset(data_vars=dict(tmax_mean=(["lat", "lon"], mean_tmax),
                                       ),
                               coords=dict(                           
                                       lat=(["lat"], ds_tmax.lat.values),
                                       lon=(["lon"], ds_tmax.lon.values),
                                       )    
                               )
        ds_output.attrs = ds_tmax.attrs
        
        Path(os.path.dirname(new_filename)).mkdir(parents=True, exist_ok=True)
        
        ds_output.to_netcdf(new_filename)
        ds_output.close()


def compute_difference(path_future_output):
    file_list = [f for f in iglob(path_future_output, recursive=True) if os.path.isfile(f)]
    
    # nr_runs = len(file_list)
    
    # ds_mean_fut = xr.open_dataset(file_list[0])
    # size = ds_mean_fut.tmax_mean.values.shape
    
    # multi_model_diff = np.empty((nr_runs, size[0], size[1]))
    
    for idx_file, fut_file in enumerate(file_list):
        his_file0 = fut_file.replace("future", "historic")
        his_file1 = his_file0.replace("ssp245/", "")
        his_file2 = his_file1.replace("ssp370/", "")
        
        directory0, directory1, directory2, directory3, directory4, str1, str2, str3, _ , _, _  = his_file2.split('_')
        
        his_file = directory0+'_'+directory1+'_'+directory2+ '_'+directory3+'_'+directory4+'_'+str1+'_'+str2+'_'+str3+'_'+'hist_r1i1p1f1_mean.nc'
        
        ds_mean_fut = xr.open_dataset(fut_file)
        ds_mean_his = xr.open_dataset(his_file)
        
        mean_fut = ds_mean_fut.tmax_mean.values
        mean_his = ds_mean_his.tmax_mean.values
        
        difference = mean_fut-mean_his
        
        
        ds_output = xr.Dataset(data_vars=dict(diff=(["lat", "lon"], difference),
                                       ),
                               coords=dict(                           
                                       lat=(["lat"], ds_mean_fut.lat.values),
                                       lon=(["lon"], ds_mean_fut.lon.values),
                                       )    
                               )
        ds_output.attrs = ds_mean_his.attrs
        
        new_filename = fut_file.replace("future", "future-historic")
        new_filename = new_filename.replace("mean", "difference")
        Path(os.path.dirname(new_filename)).mkdir(parents=True, exist_ok=True)
        
        ds_output.to_netcdf(new_filename)
        ds_output.close()
        
        # multi_model_diff[idx_file, :, :] = 
        
    
    
    # multi_model_mean = np.mean(multi_model_diff, axis=0)
    
    # return multi_model_mean


def compute_mean_per_county(path_input):

    path_soy = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/shapefiles/input/soybean_us_counties.shp'
    shapefile_soy = gpd.read_file(path_soy)
    path_corn = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/shapefiles/input/corn_us_counties.shp'
    shapefile_corn = gpd.read_file(path_corn)
    
    crops_shapefile = pd.concat([shapefile_soy,shapefile_corn]).drop_duplicates().reset_index(drop=True)
    
    dir_list = [f for f in iglob(path_input, recursive=True) if os.path.isdir(f)][2:]
    
    for directory in dir_list:
        file_list = [f for f in iglob(directory+'/*', recursive=True) if os.path.isfile(f)]
    
        multi_models = np.zeros((crops_shapefile.shape[0], len(file_list)))
        
        for file_idx, file in enumerate(file_list):
            
            # open file containing diffence in mean per model
            ds = xr.open_dataset(file)
                
            #get corresponding future file:
            fut_file0 = file.replace("future-historic", "future")
            fut_file1 = fut_file0.replace("difference", "2015-2100")
            fut_file2= fut_file1.replace("output", "input")
                
            # get affine of nc-file with rasterio
            affine = rio.open(fut_file2).transform
                
            
            # go through all geometries and compute zonal statistics
            list_model_run =[]
            for shape in crops_shapefile.geometry:
                list_model_run.append(rstats.zonal_stats(shape, ds['diff'].values.T, affine=affine, stats="mean", all_touched=True)[0]['mean'])
            
            multi_models[:, file_idx] = np.asarray(list_model_run)
        
        _, _, _, _, _, _, _, _, _, _, ssp, season = directory.split('/')
        crops_shapefile[ssp+'_'+season] = np.nanmean(multi_models, axis=1)
    
    return crops_shapefile

# """Compute mean per time period and grid cell for historic and future time period"""
# # future
# path_future = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/future/input/**/*' 
# # '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/output'
# nr_years = 40
# create_mean(path_future, nr_years, "2015-2100")


# # historic
# path_historic = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/historic/input/**/*' 
# # '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/output'

# create_mean(path_historic, nr_years, "1961-2014")


# """Compute multi model mean"""
# path_future_output = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/future/output/**/*' 
# compute_difference(path_future_output)


"""Compute value for each county"""
path_difference = '/Users/carmenst/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/CMIP6_DATA/future-historic/output/**/*' 
crops_shapefile = compute_mean_per_county(path_difference)
crops_shapefile.to_file('counties_multi_model_mean.shp')

# with rasterio.open("/path/to/raster.tif") as src:
#     affine = src.transform
#     array = src.read(1)
# df_zonal_stats = pd.DataFrame(zonal_stats(shapefile_soy, tmax, affine=affine))

# # adding statistics back to original GeoDataFrame
# gdf2 = pd.concat([gdf, df_zonal_stats], axis=1)

