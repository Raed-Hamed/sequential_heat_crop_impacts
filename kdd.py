#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:19:10 2022

@author: carmensteinmann
"""

from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd

"""Thresholds and data paths"""
#temperatur thresholds for spring and summer temperatures
t_high_spring = 24
t_high_summer = 31

#number of days prior to the harvest day
spring_start = 90
spring_end = 60
summer_start = 50
summer_end = 20

#specify temperature and soil data path and harvest date file
path_tmax =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output/tmax180')
# path_moisture = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/soil_moisture')
path_crop_map = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/crop_maps/sacks_wheat_winter.harvest.doy.I.nc')


"""Read data files"""
#grid cells with a harvest day
ds_harvest = xr.open_dataset(path_crop_map)
harvest_end_mean = ds_harvest['harvest.mean'].values
grid_cells = np.where(~np.isnan(harvest_end_mean))

#read tmax files
files = [f.name for f in path_tmax.iterdir() if f.is_file() and not f.name.startswith('.')]
files.sort()
#get data dimensions from one tmax file
ds_test = xr.open_dataset(Path(path_tmax, files[0]))
tmax_test = ds_test.tmax.values
tmax_test2 = np.empty(np.append(ds_test.tmax.values.shape, 1))
tmax_test2[:,grid_cells[0], grid_cells[1],0] = ds_test.tmax.values[:,grid_cells[0], grid_cells[1]]
lat_tmax = ds_test.lat.values
lon_tmax = ds_test.lon.values

#load temperature data for all years
nr_years_tmax= len(files)
tmax = np.empty(np.append(ds_test.tmax.values.shape, nr_years_tmax))
time = []
for idx_file, file in enumerate(files):
    ds = xr.open_dataset(Path(path_tmax, file))
    tmax[:, grid_cells[0],grid_cells[1], idx_file] = ds.tmax.values[:,grid_cells[0], grid_cells[1]]
    time.append(str(pd.to_datetime(ds.time.values[0]).year))

#load soil moisture data
# moisture = None

"""Compute kdd for summer and spring for all grid cells with a harvest date"""
sum_kdd_spring = np.empty((nr_years_tmax, 360, 720))
sum_kdd_spring[:,:,:] = np.nan
sum_kdd_summer = np.empty((nr_years_tmax, 360, 720))
sum_kdd_summer[:,:,:] = np.nan

# nr_years_moisture = 30
# moisture_spring = np.zeros((nr_years_moisture, 360, 720))
# moisture_summer = np.zeros((nr_years_moisture, 360, 720))

for idx, _ in enumerate(grid_cells[0]):
    #grid cell
    lat = grid_cells[0][idx]
    lon = grid_cells[1][idx]
    harvest_grid_cell = harvest_end_mean[lat, lon]

    #kdd spring
    day_spring_start = int(harvest_grid_cell-spring_start)
    day_spring_end = int(harvest_grid_cell-spring_end)
    kdd_spring_gridcell = tmax[day_spring_start:day_spring_end, lat, lon, :] - t_high_spring
    kdd_spring_gridcell[kdd_spring_gridcell<=0] = 0
    #kdd_spring_gridcell[np.isnan(kdd_spring_gridcell)] = 0
    sum_kdd_spring[:, lat, lon] = np.sum(kdd_spring_gridcell, axis=0)

    #kdd summer
    day_summer_start = int(harvest_grid_cell-summer_start)
    day_summer_end = int(harvest_grid_cell-summer_end)
    kdd_summer_gridcell = tmax[day_summer_start:day_summer_end, lat, lon, :] - t_high_summer
    kdd_summer_gridcell[kdd_summer_gridcell<=0] = 0
    #kdd_summer_gridcell[np.isnan(kdd_summer_gridcell)] = 0
    sum_kdd_summer[:, lat, lon] = np.sum(kdd_summer_gridcell, axis=0)
    
    # #soil moisture
    # moisture_spring_gridcell = moisture[day_spring_start:day_spring_end, lat, lon, :]
    # moisture_spring[:,lat, lon] = np.mean(moisture_spring_gridcell, axis=0)
    # moisture_summer_gridcell = moisture[day_summer_start:day_summer_end, lat, lon, :]
    # moisture_summer[:,lat, lon] = np.mean(moisture_summer_gridcell, axis=0)


"""Save output"""
output_path = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output')
filename_output = 'kdd.nc'
ds_output = xr.Dataset(data_vars=dict(kdd_spring=(["time", "lat", "lon"], sum_kdd_spring),
                               kdd_summer=(["time", "lat", "lon"], sum_kdd_summer),),
                       coords=dict(
                           time=(["time"], time),
                           lat=(["lat"], lat_tmax),
                           lon=(["lon"], lon_tmax),
                           )
                       )

ds_output.to_netcdf(Path(output_path, filename_output))
ds_output.close()


#test 
# test_kdd = sum_kdd_summer[6,107,491]
# year = 6
# lat = 107
# lon = 491
# harvest_grid_cell = harvest_end_mean[lat, lon]
# day_summer_start = int(harvest_grid_cell-summer_start)
# day_summer_end = int(harvest_grid_cell-summer_end)
# kdd_summer_gridcell = tmax[day_summer_start:day_summer_end, lat, lon, :] - t_high_summer
# kdd_summer_gridcell[kdd_summer_gridcell<=0] = 0
# kdd_summer_gridcell[np.isnan(kdd_summer_gridcell)] = 0
# size_kdd = kdd_summer_gridcell.shape
# total_kdd = np.sum(kdd_summer_gridcell, axis=0)




# ds_test_output= xr.open_dataset(Path(output_path, filename_output))



# #compute temperatures over t_high_spring and t_high_summer for one file
# kdd_spring_idx = np.where(~np.isnan(tmax_test)) and np.where((tmax_test >= t_high_spring))
# tmax_over_t_high_spring = np.empty(tmax_test.shape)
# tmax_over_t_high_spring[kdd_spring_idx] = tmax_test[kdd_spring_idx] - t_high_spring

# kdd_summer_idx = np.where(~np.isnan(tmax_test)) and np.where((tmax_test >= t_high_summer))
# tmax_over_t_high_summer = np.empty(tmax_test.shape)
# tmax_over_t_high_summer[kdd_summer_idx] = tmax_test[kdd_summer_idx] - t_high_summer

# #compute temperatures over t_high_spring and t_high_summer
# kdd_spring_idx = np.where(~np.isnan(tmax[:, grid_cells[0],grid_cells[1], :])) and np.where((tmax >= t_high_spring))
# tmax_over_t_high_spring = np.empty(tmax.shape)
# tmax_over_t_high_spring[kdd_spring_idx] = tmax[kdd_spring_idx] - t_high_spring

# kdd_summer_idx = np.where(~np.isnan(tmax)) and np.where((tmax >= t_high_summer))
# tmax_over_t_high_summer = np.empty(tmax.shape)
# tmax_over_t_high_summer[kdd_summer_idx] = tmax[kdd_summer_idx] - t_high_summer





    # kdd_spring_idx = np.where(~np.isnan(tmax[
    #     day_spring_start:day_spring_end, lat, lon, :])) and np.where((tmax[
    #         day_spring_start:day_spring_end, lat, lon, :] >= t_high_spring))
    # sum_kdd_spring[lat, lon, :] = np.sum(tmax_over_t_high_spring[
    #     int(harvest_grid_cell-spring_start):int(harvest_grid_cell-spring_end), lat, lon])
    # sum_kdd_summer[lat, lon, :] = np.sum(tmax_over_t_high_summer[
    #     int(harvest_grid_cell-summer_start):int(harvest_grid_cell-summer_end), lat, lon])


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