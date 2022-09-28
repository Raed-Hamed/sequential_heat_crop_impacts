#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:19:10 2022

@author: carmensteinmann
"""

import xarray as xr
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

#to do: differentiate start/end per grid cell
# start_spring = 100
# end_spring = 150
# start_summer = 220
# end_summer = 250
t_high_spring = 24
t_high_summer = 31


path =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/cpc/tmax')
files = [f.name for f in path.iterdir() if f.is_file() and not f.name.startswith('.')]

ds_test = xr.open_dataset(Path(path, files[0]))
tmax_test = ds_test.tmax.values
lat_test = ds_test.lat.values
lon_test = ds_test.lon.values

#for all years
# tmax = np.empty(np.append(ds_test.tmax.values.shape, len(files)))
# for idx_file, file in enumerate(files):
#     ds = xr.open_dataset(Path(path, file))
#     tmax[:,:,:,idx_file] = ds.tmax.values

kdd_spring_idx = np.where(~np.isnan(tmax_test)) and np.where((tmax_test >= t_high_spring))
#np.nonzero(seasons) and 
tmax_over_t_high_spring = np.empty(tmax_test.shape)
tmax_over_t_high_spring[kdd_spring_idx] = tmax_test[kdd_spring_idx] - t_high_spring


tmax_over_t_high_summer = np.empty(tmax_test.shape)
tmax_over_t_high_spring[kdd_spring_idx] = tmax_test[kdd_spring_idx] - t_high_spring


path_crop_maps = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/crop_maps')
files_crop_maps = [f.name for f in path_crop_maps.iterdir() if f.is_file() and not f.name.startswith('.')]
ds_harvest = xr.open_dataset(Path(path_crop_maps, 'sacks_wheat_winter.harvest.doy.nc'))
harvest_end_mean = ds_harvest['harvest.mean'].values
#harvest_end_mean[~np.isnan(harvest_end_mean)]

mask_summer = np.empty(360,720,30)
mask_summer_lat = np.ones(lat_test)
mask_summer_lon = np.ones(lon_test)

harvest_end = np.ones((360, 720), dtype=int)
harvest_end[:,:] = int(222)

crop_maps = np.ones((360, 720, 2), dtype=int)
crop_maps[:,:,0] = int(244)
crop_maps[:,:,1] = int(303)
seasons = np.empty(tmax_test.shape)
#seasons[crop_maps[:,:,0]] = 1
seasons[crop_maps] = 1
seasons[90:120, :, :] = 1


start_spring = 100
end_spring = 150
start_summer = 220
end_summer = 250


grid_cells = np.where(~np.isnan(harvest_end_mean)) 

sum_kdd_spring = np.zeros((360, 720))
spring_start = 90
spring_end = 60
summer_start = 50
summer_end = 20

for idx, _ in enumerate(grid_cells[0]):
    lat = grid_cells[0][idx]
    lon = grid_cells[1][idx]
    harvest_grid_cell = harvest_end[lat, lon]
    sum_kdd_spring[lat, lon] = np.sum(tmax_over_t_high_spring[(harvest_grid_cell-90):(harvest_grid_cell-60), lat, lon])

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