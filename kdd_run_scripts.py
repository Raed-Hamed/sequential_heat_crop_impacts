#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:06:41 2022

@author: carmensteinmann
"""
from pathlib import Path
import xarray as xr
import numpy as np
from kdd import read_t_data, kdd_gdd_per_gridcell, save_outputs

"""Thresholds and data paths"""
#temperatur thresholds for spring kdd and gdd
t_base_spring = 0
t_optimum_spring = 22
t_high_spring = 24
gdd_max_spring = t_optimum_spring-t_base_spring
thr_spring = [t_base_spring, t_optimum_spring, t_high_spring, gdd_max_spring]

#temperatur thresholds for summer kdd and gdd
t_base_summer = 10
t_optimum_summer = 21
t_high_summer = 31
gdd_max_summer = t_optimum_summer-t_base_summer
thr_summer = [t_base_summer, t_optimum_summer, t_high_summer, gdd_max_summer]

#number of days prior to the harvest day defining spring and summer
spring_start = 90
spring_end = 60
spring_dates = [spring_start, spring_end]

summer_start = 50
summer_end = 20
summer_dates = [summer_start, summer_end]

#specify temperature and soil data path and harvest date file
path_tmax =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output/tmax180')
path_tmin =  Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output/tmin180') 
path_crop_map = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/crop_maps/sacks_wheat_winter.harvest.doy.I.nc')
# path_moisture = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/soil_moisture')

path_output = Path('/Users/carmensteinmann/Documents/CLIMADA/own_projects/sequential_heat_crop_impacts/data/output')


"""Executing functions"""
#Read harvest dates and temperature data
ds_harvest = xr.open_dataset(path_crop_map)
harvest_end_mean = ds_harvest['harvest.mean'].values
grid_cells = np.where(~np.isnan(harvest_end_mean)) #grid cells with a harvest day
tmax, lat_tmax, lon_tmax, time_tmax = read_t_data(path_tmax, 'tmax', grid_cells)
tmin, lat_tmin, lon_tmin, time_tmin = read_t_data(path_tmin, 'tmin', grid_cells)

#check whether lat, lon and time is the same for tmax and tmin

[kdd_spring, kdd_summer, gdd_spring, gdd_summer] = kdd_gdd_per_gridcell(grid_cells, harvest_end_mean, 
                                                                        tmax, tmin, thr_spring, thr_summer,
                                                                        spring_dates, summer_dates)
                                                     
filename_output = 'kdd_gdd.nc'                                                     
save_outputs(path_output, filename_output, kdd_spring, kdd_summer, gdd_spring, gdd_summer, 
                  time_tmax, lat_tmax, lon_tmax) 
                                            