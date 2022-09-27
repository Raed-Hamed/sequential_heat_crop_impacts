#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:19:10 2022

@author: carmensteinmann
"""

import xarray as xr
from pathlib import Path
import numpy as np

#to do: differentiate start/end per grid cell
start_spring = 100
end_spring = 150
start_summer = 220
end_summer = 250

t_high_spring = 24
t_high_summer = 31


path =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/cpc/tmax')
files = [f.name for f in path.iterdir() if f.is_file() and not f.name.startswith('.')]

ds_test = xr.open_dataset(Path(path, files[0]))
tmax_test = ds_test.tmax.values


tmax = np.empty(np.append(ds_test.tmax.values.shape, len(files)))
for idx_file, file in enumerate(files):
    ds = xr.open_dataset(Path(path, file))
    tmax[:,:,:,idx_file] = ds.tmax.values


tmax_spring = tmax_test[start_spring:end_spring+1, :, :]
tmax_summer = tmax_test[start_summer:end_summer+1, :, :]

#per year
kdd_spring = np.where(~np.isnan(tmax_spring)) and np.where((tmax_spring >= t_high_spring)) 

tmax_over_t_high = np.empty(tmax_spring.shape)
tmax_over_t_high[kdd_spring] = tmax_spring[kdd_spring] - t_high_spring

sum_kdd_spring = np.sum(tmax_over_t_high, axis=0)

# fig = plt.figure(figsize=(10, 8))
# m = Basemap(projection='lcc', resolution='c',
#             width=8E6, height=8E6, 
#             lat_0=45, lon_0=-100,)
# m.shadedrelief(scale=0.5)
# m.pcolormesh(lon, lat, temp_anomaly,
#              latlon=True, cmap='RdBu_r')
# plt.clim(-8, 8)
# m.drawcoastlines(color='lightgray')