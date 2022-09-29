#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:19:10 2022

@author: carmensteinmann
"""

import xarray as xr
from pathlib import Path
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt

t_high_spring = 24
t_high_summer = 31

#number of days prior to the harvest day
spring_start = 90
spring_end = 60
summer_start = 50
summer_end = 20

#specify temperature data path and harvest date file
path_tmax =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/cpc/tmax_365')
path_crop_maps = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/crop_maps')
files_crop_maps = [f.name for f in path_crop_maps.iterdir() if f.is_file() and not f.name.startswith('.')]
ds_harvest = xr.open_dataset(Path(path_crop_maps, 'sacks_wheat_winter.harvest.doy.nc'))
harvest_end_mean = ds_harvest['harvest.mean'].values
grid_cells = np.where(~np.isnan(harvest_end_mean)) 

#read tmax files
files = [f.name for f in path_tmax.iterdir() if f.is_file() and not f.name.startswith('.')]
files.sort()

#get data dimensions
ds_test = xr.open_dataset(Path(path_tmax, files[0]))
tmax_test = ds_test.tmax.values
tmax_test2 = np.empty(np.append(ds_test.tmax.values.shape, 1))
tmax_test2[:,grid_cells[0], grid_cells[1],0] = ds_test.tmax.values[:,grid_cells[0], grid_cells[1]]
lat_tmax = ds_test.lat.values
lon_tmax = ds_test.lon.values

#load data for all years
nr_years= len(files)
tmax = np.empty(np.append(ds_test.tmax.values.shape, nr_years))
time = []
for idx_file, file in enumerate(files):
    ds = xr.open_dataset(Path(path_tmax, file))
    #tmax[:,:,:,idx_file] = ds.tmax.values
    tmax[:, grid_cells[0],grid_cells[1], idx_file] = ds.tmax.values[:,grid_cells[0], grid_cells[1]]
    time.append(str(pd.to_datetime(ds.time.values[0]).year))


# #compute temperatures over t_high_spring and t_high_summer for one file
# kdd_spring_idx = np.where(~np.isnan(tmax_test)) and np.where((tmax_test >= t_high_spring))
# tmax_over_t_high_spring = np.empty(tmax_test.shape)
# tmax_over_t_high_spring[kdd_spring_idx] = tmax_test[kdd_spring_idx] - t_high_spring

# kdd_summer_idx = np.where(~np.isnan(tmax_test)) and np.where((tmax_test >= t_high_summer))
# tmax_over_t_high_summer = np.empty(tmax_test.shape)
# tmax_over_t_high_summer[kdd_summer_idx] = tmax_test[kdd_summer_idx] - t_high_summer

#compute temperatures over t_high_spring and t_high_summer 
kdd_spring_idx = np.where(~np.isnan(tmax[:, grid_cells[0],grid_cells[1], :])) and np.where((tmax >= t_high_spring))
tmax_over_t_high_spring = np.empty(tmax.shape)
tmax_over_t_high_spring[kdd_spring_idx] = tmax[kdd_spring_idx] - t_high_spring

kdd_summer_idx = np.where(~np.isnan(tmax)) and np.where((tmax >= t_high_summer))
tmax_over_t_high_summer = np.empty(tmax.shape)
tmax_over_t_high_summer[kdd_summer_idx] = tmax[kdd_summer_idx] - t_high_summer


sum_kdd_spring = np.zeros((nr_years, 360, 720))
sum_kdd_summer = np.zeros((nr_years, 360, 720)) 
for idx, _ in enumerate(grid_cells[0]):
    lat = grid_cells[0][idx]
    lon = grid_cells[1][idx]
    harvest_grid_cell = harvest_end_mean[lat, lon]
    day_spring_start = int(harvest_grid_cell-spring_start)
    day_spring_end = int(harvest_grid_cell-spring_end)
    
    kdd_spring_gridcell = tmax[day_spring_start:day_spring_end, lat, lon, :] - t_high_spring
    kdd_spring_gridcell[kdd_spring_gridcell<=0] = 0
    kdd_spring_gridcell[np.isnan(kdd_spring_gridcell)] = 0
    
    sum_kdd_spring[:, lat, lon] = np.sum(kdd_spring_gridcell, axis=0)
    
    day_summer_start = int(harvest_grid_cell-summer_start)
    day_summer_end = int(harvest_grid_cell-summer_end)
    
    kdd_summer_gridcell = tmax[day_summer_start:day_summer_end, lat, lon, :] - t_high_summer
    kdd_summer_gridcell[kdd_summer_gridcell<=0] = 0
    kdd_summer_gridcell[np.isnan(kdd_summer_gridcell)] = 0
    
    sum_kdd_summer[:, lat, lon] = np.sum(kdd_summer_gridcell, axis=0)
    
    
    
    # kdd_spring_idx = np.where(~np.isnan(tmax[
    #     day_spring_start:day_spring_end, lat, lon, :])) and np.where((tmax[
    #         day_spring_start:day_spring_end, lat, lon, :] >= t_high_spring))
    # sum_kdd_spring[lat, lon, :] = np.sum(tmax_over_t_high_spring[
    #     int(harvest_grid_cell-spring_start):int(harvest_grid_cell-spring_end), lat, lon])
    # sum_kdd_summer[lat, lon, :] = np.sum(tmax_over_t_high_summer[
    #     int(harvest_grid_cell-summer_start):int(harvest_grid_cell-summer_end), lat, lon])




output_path = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output')
filename_output = 'kdd.nc'
ds_output = xr.Dataset(data_vars=dict(kdd_spring=(["time", "lat", "lon"], sum_kdd_spring), 
                               kdd_summer=(["time", "lat", "lon"], sum_kdd_summer),),
                coords=dict(
                    time=(["time"], time),
                    lat=(["lat"], lat_tmax),
                    lon=(["lon"], lon_tmax),
                    ))
                
ds_output.to_netcdf(Path(output_path, filename_output))
ds_output.close()


# ds_test_output= xr.open_dataset(Path(output_path, filename_output))


#harvest_end_mean[~np.isnan(harvest_end_mean)]

# mask_summer = np.empty(360,720,30)
# mask_summer_lat = np.ones(lat_test)
# mask_summer_lon = np.ones(lon_test)

# harvest_end = np.ones((360, 720), dtype=int)
# harvest_end[:,:] = int(222)

# crop_maps = np.ones((360, 720, 2), dtype=int)
# crop_maps[:,:,0] = int(244)
# crop_maps[:,:,1] = int(303)
# seasons = np.empty(tmax_test.shape)
# #seasons[crop_maps[:,:,0]] = 1
# seasons[crop_maps] = 1
# seasons[90:120, :, :] = 1






# tmax_spring = tmax_test[start_spring:end_spring+1, :, :]
# tmax_summer = tmax_test[start_summer:end_summer+1, :, :]

# harvest_day = 300
# summer_duration = 30
# tmax_spring = tmax_test[harvest_day-summer_duration:harvest_day, :, :]
# tmax_summer = tmax_test[start_summer:end_summer+1, :, :]

#per year
# kdd_spring_idx = np.where(~np.isnan(tmax_spring)) and np.where((tmax_spring >= t_high_spring)) 



# sum_kdd_spring = np.sum(tmax_over_t_high, axis=0)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# m = Basemap(projection='lcc', resolution='c',
#             width=8E6, height=8E6, 
#             lat_0=45, lon_0=-100,)
# m.shadedrelief(scale=0.5)
# m.pcolormesh(lon, lat, temp_anomaly,
#               latlon=True, cmap='RdBu_r')
# plt.clim(-8, 8)
# m.drawcoastlines(color='lightgray')