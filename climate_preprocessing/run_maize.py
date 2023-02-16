#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:06:41 2022

@author: carmensteinmann
"""
from pathlib import Path
import xarray as xr
import numpy as np
from functions import read_t_data, read_moisture, from_planting_date, kdd_gdd_per_months, save_outputs, from_harvest_date

"""Thresholds and data paths"""
#temperatur thresholds for spring kdd and gdd
t_base_spring = 10
t_optimum_spring = 29
t_high_spring = 29
gdd_max_spring = t_optimum_spring-t_base_spring
thr_spring = [t_base_spring, t_optimum_spring, t_high_spring, gdd_max_spring]
        
#temperatur thresholds for summer kdd and gdd
t_base_summer = 10
t_optimum_summer = 29
t_high_summer = 29
gdd_max_summer = t_optimum_summer-t_base_summer
thr_summer = [t_base_summer, t_optimum_summer, t_high_summer, gdd_max_summer]

#number of days after planting date defining spring and summer
spring_start = 0
spring_end = 60
summer_start = 60
summer_end = 120

# #number of days after planting date defining spring and summer
# harvest_day = 160
# spring_start = harvest_day-20
# spring_end = harvest_day-50
# summer_start = harvest_day-50
# summer_end = harvest_day-80
# filename_output = 'maize_30days.nc'  

# harvest_day = 160
# spring_start = harvest_day-0
# spring_end = harvest_day-60
# summer_start = harvest_day-60
# summer_end = harvest_day-120
# filename_output = 'maize_60days.nc' 


spring_dates = [spring_start, spring_end]
summer_dates = [summer_start, summer_end]


#spring and summer defined by specific months
spring_dates_m = [90, 151] #April and May
# summer_dates_m = [152, 212] # June and July
# spring_dates_m = [121, 151] # May
# summer_dates_m = [182, 212] # July
summer_dates_m = [182, 243] # July and August

#specify temperature and soil data path and harvest date file
path_plant_day = Path('/Volumes/Files/WCR/2022/sequential_heat/data/crop_maps/plant_doy_mean_Maize.crop.calendar.nc')
path_harvest_day = Path('/Volumes/Files/WCR/2022/sequential_heat/data/crop_maps/harvest_doy_mean_Maize.crop.calendar.nc')

path_tmax =  Path('/Volumes/Files/WCR/2022/sequential_heat/data/output/tmax180')
path_tmin =  Path('/Volumes/Files/WCR/2022/sequential_heat/data/output/tmin180') 
path_moisture = Path('/Volumes/Files/WCR/2022/sequential_heat/data/SMroot_1980-2021_GLEAM_v3.6a_daily_remap05.nc')

path_output = Path('/Volumes/Files/WCR/2022/sequential_heat/data/output')


"""Executing functions"""
#Read harvest dates, temperature and soil moisture data
ds_planting = xr.open_dataset(path_plant_day)
regridded_planting = ds_planting.coarsen(longitude=6).mean().coarsen(latitude=6).mean()
planting = regridded_planting['plant'].values
#grid_cells = np.where(~np.isnan(planting)) #grid cells with a planting day

ds_harvest = xr.open_dataset(path_harvest_day)
regridded_harvest = ds_harvest.coarsen(longitude=6).mean().coarsen(latitude=6).mean()
harvest = regridded_harvest['harvest'].values
growing_time = harvest-planting

grid_cells = (np.where(~np.isnan(growing_time)) and np.where(growing_time>0)) #grid cells with a planting date earlier in the year than harvesting


tmax, lat_tmax, lon_tmax, time_tmax = read_t_data(path_tmax, 'tmax', grid_cells)
tmin, lat_tmin, lon_tmin, time_tmin = read_t_data(path_tmin, 'tmin', grid_cells)
moisture, lat_moisture, lon_moisture, time_moisture = read_moisture(path_moisture, grid_cells)


#check whether lat, lon and time is the same for tmax, tmin and soil moisture


# #kdd, gdd and mean sm for grid cell specific harvest days
# [kdd_spring, kdd_summer, gdd_spring, gdd_summer, 
#  mean_tmax_spring, mean_tmax_summer, moisture_spring, moisture_summer] = from_planting_date(grid_cells, planting,
#                                                                                                tmax, tmin, moisture, 
#                                                                                                thr_spring, thr_summer,
#                                                                                                spring_dates, summer_dates)
                                                                                            
# [kdd_spring, kdd_summer, 
#   gdd_spring, gdd_summer, 
#   mean_tmax_spring, mean_tmax_summer,
#   moisture_spring, moisture_summer] = from_harvest_date(grid_cells, harvest, 
#                                                           tmax, tmin, moisture, 
#                                                           thr_spring, thr_summer,
#                                                           spring_dates, summer_dates)
                                                   
# save_outputs(path_output, filename_output, kdd_spring, kdd_summer, gdd_spring, gdd_summer, 
#               mean_tmax_spring, mean_tmax_summer, moisture_spring, moisture_summer, time_tmax, lat_tmax, lon_tmax) 

#kdd, gdd and mean sm for fixed_months
[kdd_spring_m, kdd_summer_m, 
  gdd_spring_m, gdd_summer_m, 
  mean_tmax_spring_m, mean_tmax_summer_m,
  moisture_spring_m, moisture_summer_m] = kdd_gdd_per_months(grid_cells, tmax, tmin, moisture, 
                                                              thr_spring, thr_summer,
                                                              spring_dates_m, summer_dates_m)                                                          
filename_output = 'maize_AM_JA.nc'                                                     
save_outputs(path_output, filename_output, kdd_spring_m, kdd_summer_m, gdd_spring_m, gdd_summer_m, 
             mean_tmax_spring_m, mean_tmax_summer_m, moisture_spring_m, moisture_summer_m, time_tmax, lat_tmax, lon_tmax) 


"""Plotting"""
from functions import plot_map
plot_map(mean_tmax_spring_m[3,:,:], lat_tmax, lon_tmax)


# year = 0
# lat = 82
# lon = 140
# season = 'summer'
# filename = 'kdd_gdd_vis.pdf'
# plot_gdd_kdd(year, lat, lon, tmax, tmin, harvest_end_mean, thr_summer, summer_dates, season, filename)


                 